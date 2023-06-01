from magnetoget import download_magnetogram, reduce_fits_image

cr = 2219

data_dir = f"/Users/cgilbert/vscode/Fluxon-Scripts-Gilly"

params_path = data_dir + "magnetic_target.params"

(hmi_path, mdi_path) = download_magnetogram(cr=cr, date=None, data_dir=data_dir)
reduce_fits_image(hmi_path, target_resolution=None, reduction_amount=10)
# reduce_fits_image(mdi_path, target_resolution=None, reduction_amount=10)
