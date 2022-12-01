
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry, ApertureStats
from astropy.io import fits

from pathlib import Path
from dateutil.parser import parse



def do_photometry_with_timestamps(
        images,
        start_object_x,
        start_object_y,
        start_object_t,
        end_object_x,
        end_object_y,
        end_object_t,
        aperture_=15,
        annulus_inner=25,
        annulus_outer=45,
        dry_run=False
        ):
    images = list(images)
    
    deltax_total = end_object_x - start_object_x
    deltay_total = end_object_y - start_object_y
    end_object_t, start_object_t = parse(end_object_t).timestamp(), parse(start_object_t).timestamp()

    delta_t_total = end_object_t - start_object_t

    velocity_x = deltax_total / delta_t_total
    velocity_y = deltay_total / delta_t_total

    photometry_target = pd.DataFrame(columns = ["num", "id", "xcenter", "ycenter", "aperture_sum", "total_bkg", "aperture_sum_bkgsub"])
    photometry_star = pd.DataFrame(columns = ["num", "id", "xcenter", "ycenter", "aperture_sum", "total_bkg", "aperture_sum_bkgsub"])

    for i, img_path in enumerate(images):
        img_open = fits.open(Path(img_path))
        img = img_open[0].data

        t_obj = parse(img_open[0].header['DATE-OBS']).timestamp()

        delta_t_obj = t_obj - start_object_t
        object_position_x = start_object_x + velocity_x * delta_t_obj
        object_position_y = start_object_y + velocity_y * delta_t_obj
        
        #target_position = (348, 971)
        #positions = [target_position, (object_position_x, object_position_y)]
        positions = [(object_position_x, object_position_y)] # ADD OTHERS?

        aperture = CircularAperture(positions, r = aperture_)
        annulus_aperture = CircularAnnulus(positions, r_in = annulus_inner, r_out = annulus_outer)
        
        if dry_run:
            if not i % (len(images) // 5):
                fig, ax = plt.subplots()
                ax.invert_yaxis()

                ax.imshow(img, cmap='gray', vmin=0, vmax=255)
                aperture.plot(ax, color='red', lw=1, label='Photometry aperture')
                annulus_aperture.plot(ax, color='green', lw=1, label='Photo annulus')

                plt.show(block=True)
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
            photometry_target.loc[i] = phot_dataframe.iloc[0]
            # photometry_star.loc[i] = phot_dataframe.iloc[1]
        ## add more here if there are more than 2 apertures
        
    combined_data = pd.DataFrame(columns = ["target_subtracted_counts", "object_subtracted_counts", "num"])
    combined_data["target_subtracted_counts"] = photometry_target["aperture_sum_bkgsub"]
    combined_data["object_subtracted_counts"] = photometry_star["aperture_sum_bkgsub"]
    combined_data["differential"] = combined_data["object_subtracted_counts"] - combined_data["target_subtracted_counts"]
    combined_data["uncalibrated_object_mag"] = -2.5 * np.log10(combined_data["object_subtracted_counts"])
    combined_data["uncalibrated_target_mag"] = -2.5 * np.log10(combined_data["target_subtracted_counts"])
    combined_data["differential_mag"] = combined_data["uncalibrated_object_mag"] - combined_data["uncalibrated_target_mag"]
    for i in range(len(combined_data)):
        combined_data["num"].loc[i] = i
    #plt.scatter(combined_data["num"], combined_data["differential"])

    return combined_data