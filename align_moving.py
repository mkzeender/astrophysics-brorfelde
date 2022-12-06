
from collections import namedtuple
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry, ApertureStats
from astropy.io import fits

from pathlib import Path
from dateutil.parser import parse
import mplcursors



def get_coords_of_point(images, vmin=0, vmax=255):
    if isinstance(images, str):
        images = [images]

    with plt.ion():
        data = []
        for image in images:
            with fits.open(image) as img_open:
                img = img_open[0].data
                t = img_open[0].header['DATE-OBS']

                data.append(ApPhotIm(t, img,)) #CircularAperture([], 13*3)))

                # def print_coords(sel):
                #     x, y = sel.target
                #     print(repr(TargetLoc(x, y, t)))
                #     sel.annotation.set_text(f'{sel.target[0]}, {sel.target[1]}, {repr(t)}')

                # fig, ax = plt.subplots()

                # _im = ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, origin='lower')

                # mplcursors.cursor(_im).connect("add", print_coords)

                # plt.show(block=True)
        
        MotionFigure(data=data, vmin=vmin, vmax=vmax).show()

            

class DryRunDone(BaseException):
    pass

@dataclass(order=True)
class ApPhotIm:
    t: datetime
    img: np.ndarray
    ap: Optional[CircularAperture] = None
    ann: Optional[CircularAnnulus] = None

@dataclass
class TargetLoc:
    x: float
    y: float
    t: Optional[str] = None

    @property
    def t_datetime(self):
        return parse(self.t)
    
    def timestamp(self):
        return self.t_datetime.timestamp()

class MotionFigure:
    def __init__(self, data:List[ApPhotIm]=[], vmin=0, vmax=255) -> None:
        self.vmin = vmin
        self.vmax = vmax

        self.fig: Figure = plt.figure()

        self.ax: plt.Axes = self.fig.add_axes([0, 0, 1, 1])

        self.data = data.copy()
        self.data.sort()

        self.i = 0

        self.fig.canvas.mpl_connect('key_press_event', self._on_key)

        self.curs: Optional[mplcursors.Cursor] = None

    def _on_key(self, event):
        if event.key == 'n':
            self.next()
        elif event.key == 'b':
            self.prev()

        
    
    def next(self):
        self.i += 1
        self.i %= len(self.data)
        self._show()

    def prev(self):
        self.i -= 1
        self.i %= len(self.data)
        self._show()

    def _show(self):
        dat = self.data[self.i]

        self.ax.remove()
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.imshow(dat.img, cmap='gray', vmax=self.vmax, vmin=self.vmin)
        if dat.ap is not None:
            dat.ap.plot(self.ax, color='red')
        if dat.ann is not None:
            dat.ann.plot(self.ax, color='green')

        if self.curs is not None:
            self.curs.remove()

        self.curs = mplcursors.cursor(self.ax)
        self.curs.connect("add", self._print_coords)

        self.fig.canvas.draw_idle()

    def _print_coords(self, sel):
        x, y = sel.target
        print(repr(TargetLoc(x, y, self.data[self.i].t)))
        sel.annotation.set_text(f'{sel.target[0]}, {sel.target[1]}, {self.data[self.i].t}')

    def show(self):
        self._show()
        plt.show(block=True)

        

def do_photometry_with_timestamps(
        images: List[str],
        start_object: TargetLoc,
        end_object: TargetLoc,
        aperture_=15,
        annulus_inner=25,
        annulus_outer=45,
        dry_run=False,

        reference_positions: List[TargetLoc] = None,

        vmin = 0,
        vmax = 255
        ):
    images = list(images)
    
    deltax_total = end_object.x - start_object.x
    deltay_total = end_object.y - start_object.y
    end_object_t, start_object_t = end_object.timestamp(), start_object.timestamp()

    delta_t_total = end_object_t - start_object_t

    velocity_x = deltax_total / delta_t_total
    velocity_y = deltay_total / delta_t_total

    photometry_target = pd.DataFrame(columns = ["num", "id", "xcenter", "ycenter", "aperture_sum", "total_bkg", "aperture_sum_bkgsub", 'timestamp'])
    photometry_star_list = [
        pd.DataFrame(columns = ["num", "id", "xcenter", "ycenter", "aperture_sum", "total_bkg", "aperture_sum_bkgsub"])
        for _ in reference_positions
    ]

    
    
    _dry_run_info = []

    for i, img_path in enumerate(images):
        with fits.open(Path(img_path)) as img_open:
            img = img_open[0].data

            t_obj_readable = parse(img_open[0].header['DATE-OBS'])
            t_obj = t_obj_readable.timestamp()

        delta_t_obj = t_obj - start_object_t
        object_position_x = start_object.x + velocity_x * delta_t_obj
        object_position_y = start_object.y + velocity_y * delta_t_obj
        
        #target_position = (348, 971)
        positions = [(p.x, p.y) for p in reference_positions] + [(object_position_x, object_position_y)]
        # positions = [(object_position_x, object_position_y)] # ADD OTHERS?

        aperture = CircularAperture(positions, r = aperture_)
        annulus_aperture = CircularAnnulus(positions, r_in = annulus_inner, r_out = annulus_outer)
        
        if dry_run:
            # if not i % (len(images) // 5):
            #     fig, ax = plt.subplots()


            #     ax.imshow(img, cmap='gray', vmin=0, vmax=255, origin='upper')
            #     aperture.plot(ax, color='red', lw=1, label='Photometry aperture')
            #     annulus_aperture.plot(ax, color='green', lw=1, label='Photo annulus')

            #     plt.show(block=True)
            _dry_run_info.append(ApPhotIm(t=t_obj_readable, img=img, ap=aperture, ann=annulus_aperture))

        else:
            phot_table = aperture_photometry(img, aperture)
            aperstats = ApertureStats(img, annulus_aperture)
            bkg_mean = aperstats.mean
            aperture_area = aperture.area_overlap(img)
            total_bkg = bkg_mean * aperture_area
            phot_table['total_bkg'] = total_bkg
            phot_bkgsub = phot_table['aperture_sum'] - total_bkg
            phot_table['aperture_sum_bkgsub'] = phot_bkgsub

            for col in phot_table.colnames:
                phot_table[col].info.format = '%.8g'  # for consistent table output
            

            phot_dataframe = pd.DataFrame(np.array(phot_table))
            for j in range(len(reference_positions)):
                photometry_star_list[j].loc[i] = phot_dataframe.iloc[j]
            photometry_target.loc[i] = phot_dataframe.iloc[j+1]
            photometry_target.loc[i, 'timestamp'] = t_obj_readable
        ## add more here if there are more than 2 apertures

    #photometry_star_mean = pd.DataFrame(columns = ["num", "id", "xcenter", "ycenter", "aperture_sum", "total_bkg", "aperture_sum_bkgsub"])

    if dry_run:
        motion_fig = MotionFigure(data=_dry_run_info, vmin=vmin, vmax=vmax)
        motion_fig.show()

        raise DryRunDone('#not_an_error')

    photometry_star_mean = photometry_star_list[0]
    # calculate the mean values for all the stars
    for data_frame in photometry_star_list[1:]:
        photometry_star_mean += data_frame
    
    
    photometry_star_mean /= len(photometry_star_list)

    combined_data = pd.DataFrame(columns = ["target_subtracted_counts", "star_subtracted_counts", "num", 'timestamp'])
    combined_data["target_subtracted_counts"] = photometry_target["aperture_sum_bkgsub"]
    combined_data["star_subtracted_counts"] = photometry_star_mean["aperture_sum_bkgsub"]
    combined_data["differential"] = combined_data["target_subtracted_counts"] - combined_data["star_subtracted_counts"]
    combined_data["uncalibrated_star_mag"] = -2.5 * np.log10(combined_data["star_subtracted_counts"])
    combined_data["uncalibrated_target_mag"] = -2.5 * np.log10(combined_data["target_subtracted_counts"])
    combined_data["differential_mag"] = combined_data["uncalibrated_target_mag"] - combined_data["uncalibrated_star_mag"]

    combined_data['timestamp'] = photometry_target['timestamp']
    for i in range(len(combined_data)):
        combined_data.loc[i, "num"] = i
    #plt.scatter(combined_data["num"], combined_data["differential"])

    return combined_data