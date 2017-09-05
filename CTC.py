from astropy.io import ascii

def main():
	array_limit_file = 'CHARA_Global_Limits.txt'
	combiner_limit_file = 'Beam_Combiner_Limits.txt'
	dec_lim,Rmag_lim = get_array_limits(array_limit_file)
	combiner_dict = get_combiner_limits(combiner_limit_file)
	
	for i in combiner_dict:
		print i,combiner_dict[i]
	
def get_array_limits(array_limit_file):
	"""
	This function reads the text file that has the Array's global limits (Declination and R-band for Tracking)
	Inputs:
		array_limit_file - The filename for the Array's global limits. Must be in the same directory as CTC.py
	Outputs:
		dec_lim - The Array's declination limit
		Rmag_lim - The Array's limiting R-band magnitude for tiptilt
	"""
	data = ascii.read(array_limit_file)			#Reads the file
	params = list(data['Parameter'])			#Makes a list of the parameter names
	limits = list(data['Limit'])				#Makes a list of the parameter values
	for i in range(len(params)):				#Goes through the params
		if params[i] == 'Declination (deg)':		#If the current param is Declination
			dec_lim = float(limits[i])					#Then set dec_lim to the current value
		if params[i] == 'Tracking R-band':			#If the current param is R-band
			Rmag_lim = float(limits[i])					#Then set Rmag_lim to the current value
	return dec_lim,Rmag_lim
def get_combiner_limits(combiner_limit_file):
	"""
	This function reads the text file that has the Array's global limits (Declination and R-band for Tracking)
	Inputs:
		combiner_limit_file - The filename for the Beam Combiners' limits. Must be in the same directory as CTC.py
	Outputs:
		combiner_dict - A dictionary with keys being each combiner/bandpass combo and is associated with a list
			of values : [Bandpass, Mag_Limit, Diam_Lim_mas, Img_Lim_mas] - The bandpass, magnitude faint limit,
			limit for measuring diameters, and limit for doing imaging.
	"""
	data = ascii.read(combiner_limit_file)		#Reads the file
	combiners = list(data['Combiner'])
	bandpasses = list(data['Bandpass'])
	mag_lims = list(data['Mag_Limit'])
	diam_lims = list(data['Diam_Lim_mas'])
	img_lims = list(data['Img_Lim_mas'])
	combiner_dict = dict()
	for i in range(len(combiners)):
		this_bandpass = bandpasses[i]
		if mag_lims[i] == 'N/A':
			this_mag_lim = ''
		else:
			this_mag_lim = float(mag_lims[i])
		if diam_lims[i] == 'N/A':
			this_diam_lim = ''
		else:
			this_diam_lim = float(diam_lims[i])
		if img_lims[i] == 'N/A':
			this_img_lim = ''
		else:
			this_img_lim = float(img_lims[i])
		combiner_dict[combiners[i]] = [this_bandpass,this_mag_lim,this_diam_lim,this_img_lim]
	return combiner_dict

if __name__=="__main__":
	main()