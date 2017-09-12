from astropy.io import ascii
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import sys

def main():
	#dec = -20.
	#Rmag = 11.
	#Hmag = 8.
	#Kmag = 8.
	#diam = 0.6
	
	#star = 'kap And'
	star=''
	for i in range(len(sys.argv)):
		if i != 0:
			star+=sys.argv[i]
			if i != len(sys.argv)-1:
				star+=' '
	
	v = Vizier(columns=['*','_r'])
	table = v.query_object(star)['II/346/jsdc_v2']
	#print table
	
	rs = list(table['_r'])

	i = np.where(rs == min(rs))[0][0]
	
	coord = table['RAJ2000'][i]+' '+table['DEJ2000'][i]
	coord = SkyCoord(coord, unit=(u.hourangle,u.deg))
	dec  = coord.dec.degree
	Rmag = table['Rmag'][i]
	Hmag = table['Hmag'][i]
	Kmag = table['Kmag'][i]
	diam = table['LDD'][i]
	
	
	print 'Input Target Name: {}'.format(star)
	print 'Found Target Name: {}'.format(table['Name'][i])
	print 'Dec: {} deg, Rmag: {} mag, Hmag: {} mag, Kmag: {} mag, diam: {} mas'.format(dec,Rmag,Hmag,Kmag,diam) 
	
	
	do_classic = True
	do_climb = False
	do_jouflu = False
	do_mirc = False
	do_pavo = False
	do_vega = False
	
	array_limit_file = 'CHARA_Global_Limits.txt'
	combiner_limit_file = 'Beam_Combiner_Limits.txt'
	dec_lim,Rmag_lim = get_array_limits(array_limit_file)
	dec_flag,track_flag = check_target_global(dec,Rmag,dec_lim,Rmag_lim)
	
	print '=========='
	if 'R' not in [dec_flag,track_flag]:
		if dec_flag == 'Y' : print 'Target is low. It may be difficult to get good data.'
		if track_flag == 'Y' : print 'Target is very faint to track on. It may be difficult to get good data.'
		if 'Y' not in [dec_flag,track_flag] : print 'Target can be tracked!'
		combiner_dict = get_combiner_limits(combiner_limit_file)
		print '=========='
		if do_classic:
			classic_Hmag_flag,classic_Hmag_diam_flag,classic_Kmag_flag,classic_Kmag_diam_flag = check_target_classic(diam,combiner_dict,Hmag,Kmag)
			if classic_Hmag_flag == 'Y' : print 'Target is faint for Classic using H-band. It may be difficult to get good data.'
			if classic_Hmag_diam_flag == 'Y' : print 'Target is small for Classic using H-band. It may be difficult to get good data.'
			if classic_Hmag_flag == 'R' : print 'Target too faint for Classic using H-band.'
			if classic_Hmag_diam_flag == 'R' : print 'Target too small for Classic using H-band.'
			if 'Y' not in [classic_Hmag_flag,classic_Hmag_diam_flag] and 'R' not in [classic_Hmag_flag,classic_Hmag_diam_flag] : print 'Target can be resolved with Classic using H-band!'
			print '=========='
			if classic_Kmag_flag == 'Y' : print 'Target is faint for Classic using K-band. It may be difficult to get good data.'
			if classic_Kmag_diam_flag == 'Y' : print 'Target is small for Classic using K-band. It may be difficult to get good data.'
			if classic_Kmag_flag == 'R' : print 'Target too faint for Classic using K-band.'
			if classic_Kmag_diam_flag == 'R' : print 'Target too small for Classic using K-band.'
			if 'Y' not in [classic_Kmag_flag,classic_Kmag_diam_flag] and 'R' not in [classic_Kmag_flag,classic_Kmag_diam_flag] : print 'Target can be resolved with Classic using K-band!'
	else:
		if dec_flag == 'R' : print 'Target is too low.'
		if track_flag == 'R' : print 'Target is too faint to track on.'
		if dec_flag == 'Y' : print 'Target is low. It may be difficult to get good data.'
		if track_flag == 'Y' : print 'Target is very faint to track on. It may be difficult to get good data.'
		
	
	
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
	typ_lim = list(data['Typical Limit'])		#Makes a list of the parameter values
	best_lim = list(data['Best Limit'])			#Makes a list of the parameter values
	for i in range(len(params)):				#Goes through the params
		if params[i] == 'Declination (deg)':				#If the current param is Declination
			dec_lim = [float(typ_lim[i]),float(best_lim[i])]	#Then set dec_lim to the current value
		if params[i] == 'Tracking R-band':					#If the current param is R-band
			Rmag_lim = [float(typ_lim[i]),float(best_lim[i])]	#Then set Rmag_lim to the current value
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
	typ_mag_lims = list(data['Typ_Mag_Lim'])
	best_mag_lims = list(data['Best_Mag_Lim'])
	diam_lims = list(data['Diam_Lim_mas'])
	img_lims = list(data['Img_Lim_mas'])
	combiner_dict = dict()
	for i in range(len(combiners)):
		this_bandpass = bandpasses[i]
		this_typ_mag_lim = float(typ_mag_lims[i])
		this_best_mag_lim = float(best_mag_lims[i])
		if diam_lims[i] == 'N/A':
			this_diam_lim = ''
		else:
			this_diam_lim = float(diam_lims[i])
		if img_lims[i] == 'N/A':
			this_img_lim = ''
		else:
			this_img_lim = float(img_lims[i])
		combiner_dict[combiners[i]] = [this_bandpass,this_typ_mag_lim,this_best_mag_lim,this_diam_lim,this_img_lim]
	return combiner_dict
def check_target_global(dec,Rmag,dec_lim,Rmag_lim):
	if dec > dec_lim[0]:
		dec_flag = 'G'
	elif dec > dec_lim[1]:
		dec_flag = 'Y'
	else:
		dec_flag = 'R'

	if Rmag < Rmag_lim[0]:
		track_flag = 'G'
	elif Rmag < Rmag_lim[1]:
		track_flag = 'Y'
	else:
		track_flag = 'R'
	return dec_flag,track_flag

def check_target_classic(diam,combiner_dict,Hmag=100.,Kmag=100.):
	if Hmag != 100.:
		classic_Hmag_typ_lim = combiner_dict['Classic_H'][1]
		classic_Hmag_best_lim = combiner_dict['Classic_H'][2]
		
		if Hmag < classic_Hmag_typ_lim:
			classic_Hmag_flag = 'G'
		elif Hmag < classic_Hmag_best_lim:
			classic_Hmag_flag = 'Y'
		else:
			classic_Hmag_flag = 'R'
			
		classic_Hmag_diam_lim = combiner_dict['Classic_H'][3]
		
		if diam > classic_Hmag_diam_lim+0.1:
			classic_Hmag_diam_flag = 'G'
		elif diam > classic_Hmag_diam_lim:
			classic_Hmag_diam_flag = 'Y'
		else:
			classic_Hmag_diam_flag = 'R'
	else:
		classic_Hmag_flag = ''
		classic_Hmag_diam_flag = ''
		
	if Kmag != 100.:
		classic_Kmag_typ_lim = combiner_dict['Classic_K'][1]
		classic_Kmag_best_lim = combiner_dict['Classic_K'][2]
		
		if Kmag < classic_Kmag_typ_lim:
			classic_Kmag_flag = 'G'
		elif Kmag < classic_Kmag_best_lim:
			classic_Kmag_flag = 'Y'
		else:
			classic_Kmag_flag = 'R'
			
		classic_Kmag_diam_lim = combiner_dict['Classic_K'][3]
		
		if diam > classic_Kmag_diam_lim+0.1:
			classic_Kmag_diam_flag = 'G'
		elif diam > classic_Kmag_diam_lim:
			classic_Kmag_diam_flag = 'Y'
		else:
			classic_Kmag_diam_flag = 'R'
	else:
		classic_Kmag_flag = ''
		classic_Kmag_diam_flag = ''
	return classic_Hmag_flag,classic_Hmag_diam_flag,classic_Kmag_flag,classic_Kmag_diam_flag

def check_target_old(dec,Rmag,Hmag,Kmag,diam,dec_lim,Rmag_lim,combiner_dict):
	if dec > dec_lim+20.:
		dec_flag = 'G'
	elif dec > dec_lim:
		dec_flag = 'Y'
	else:
		dec_flag = 'R'

	if Rmag < Rmag_lim-1.:
		track_flag = 'G'
	elif Rmag < Rmag_lim:
		track_flag = 'Y'
	else:
		track_flag = 'R'
	
	combiner_mag_flags = dict()
	combiner_diam_flags = dict()
	combiner_img_flags = dict()
	

	for i in combiner_dict:
		if combiner_dict[i][0] == 'R':
			this_mag = Rmag
		if combiner_dict[i][0] == 'H':
			this_mag = Hmag
		if combiner_dict[i][0] == 'K':
			this_mag = Kmag
		
		this_mag_lim = combiner_dict[i][1]
		
		if this_mag < this_mag_lim-1.:
			combiner_mag_flags[i] = 'G'
		elif this_mag < this_mag_lim:
			combiner_mag_flags[i] = 'Y'
		else:
			combiner_mag_flags[i] = 'R'
		
		this_diam_lim = combiner_dict[i][2]
		
		if this_diam_lim != '':
			if diam > this_diam_lim-0.1:
				combiner_diam_flags[i] = 'G'
			elif diam > this_diam_lim:
				combiner_diam_flags[i] = 'Y'
			else:
				combiner_diam_flags[i] = 'R'
		else:
			combiner_diam_flags[i] = ''
				
		this_img_lim = combiner_dict[i][3]
		
		if this_img_lim != '':
			if diam > this_img_lim-0.1:
				combiner_img_flags[i] = 'G'
			elif diam > this_img_lim:
				combiner_img_flags[i] = 'Y'
			else:
				combiner_img_flags[i] = 'R'
		else:
			combiner_img_flags[i] = ''
		all_flags = [dec_flag,track_flag,combiner_mag_flags[i],combiner_diam_flags[i],combiner_img_flags[i]]
		if 'R' in all_flags:
			print 'Target cannot be observed with {}'.format(i)
		elif 'Y' in all_flags:
			print 'Target is close to limits for {}. Be wary.'.format(i)
		else:
			print 'Target is good for {}!!'.format(i)
	
if __name__=="__main__":
	main()