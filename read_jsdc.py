from astropy.io import fits

filename='jsdc2.fit'

hdulist = fits.open(filename)
for i in range(len(hdulist)):		#for all the "pages" of the fits file...
	print '***{}****'.format(i)		#print the "page" number
	print hdulist[i].header			#print the header for "page" i
	print hdulist[i].data			#print the data for "page" i

hdulist.close()
