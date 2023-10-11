import numpy as np
# import astropy as ap


def ind_split(time, gapsize=0.5):


    """
    Find indices on which data gap is larger than required.

    @time -must be in chronological order (array)
    @gapsize - minimum gap size required to split data

    """
    dataa = time[:-1]
    datab = time[1:]
    diff = datab-dataa
#    elements = np.linspace(0,len(diff)-1, len(diff)).astype(int)
#    cuts = elements[diff>gap]
    cuts = np.where(diff>gapsize)[0]
    cuts = np.append(cuts, len(time)-1)
    return cuts


def array_split(data, indices):

    """
    Split mxn array at given indices in n.
    @data - mxn data array
    @indices- indices where data should be cut

    """
#    indices = np.append(indices, int(len(data)-1))
#    cuts = np.append(cuts, int(len(data)-1))
    datasn = []

    init = 0
    for x in indices:
        datan = data[init:x+1]
        datasn.append(datan)
        init = x+1

    return datasn


def asas_clean(filename):

    # load asas file in standard format and take only best quality data

    data = np.loadtxt(filename, dtype='str')
    datan = data[(data[:,11] == 'A')]
    if len(datan) < 30:
        print(" No good quality data exists")
        datas = np.array([])
        return datas
    else:
    # pick best aperture
        col_num = len(datan[0])
        if col_num == 13:
            datan = np.delete(datan, [col_num-2,col_num-1], axis=1)
        else:
            datan = np.delete(datan, [col_num-3,col_num-2,col_num-1], axis=1)
        datan = datan.astype('float')

    #split into groups of mag and accompanying error
    #the numbers match the names in the file which is why they appear out of order

        data3 = np.column_stack((datan[:,1], datan[:,6]))
        data0 = np.column_stack((datan[:,2], datan[:,7]))
        data1 = np.column_stack((datan[:,3], datan[:,8]))
        data2 = np.column_stack((datan[:,4], datan[:,9]))
        data4 = np.column_stack((datan[:,5], datan[:,10]))

        errs = [np.median(data3[:,1]),np.median(data0[:,1]),np.median(data1[:,1]),np.median(data2[:,1]),np.median(data4[:,1])]
        datasets = [data3, data0, data1, data2, data4]
        min_ind = errs.index(min(errs))
        mag_data = datasets[min_ind]
    # recombine and sort in time

        datan = np.column_stack((datan[:,0], mag_data))
        y = np.argsort(datan[:,0])
        datan = datan[y]

    # split into chunks and remove ones that are far from the median
    #@todo make this more sophisticated. It's not necessary for the quick look
    #so right now this is good enough. numbers used here are purely empirical

        inds = ind_split(datan[:,0], 90)
        chunks = array_split(datan, inds)

        med_tot = np.median(datan[:,1][datan[:,0]>3000])

        chunks_new = []

        for x in range(len(chunks)):

            chunk = chunks[x]
            diff = med_tot - np.median(chunk[:,1])

            if abs(diff) < 1.5:

                chunks_new.append(chunk)

        datas = np.vstack(chunks_new)

# remove values with errors much larger than the mean error

        med_err = np.median(datas[:,2])
        std = np.std(datas[:,2])

        max_err = med_err+2*std

        datas = datas[datas[:,2] < max_err]

# remove large outliers

        med = np.median(datas[:,1])
        std = np.std(datas[:,1])
        max_m = med+4*std
        min_m = med-4*std
        datas = datas[(datas[:,1] < max_m) & (datas[:,1] > min_m)]

        return datas


def kws_clean(filename):

    # assumes you are using the files created from data_grabber

    data = np.loadtxt(filename, usecols=[0,1,2])

    # remove values with errors much larger than the mean error

    med_err = np.median(data[:,2])
    std = np.std(data[:,2])

    max_err = med_err+2*std

    data = data[data[:, 2] < max_err]

# remove large outliers

    med = np.median(data[:,1])
    std = np.std(data[:,1])
    max_m = med+4*std
    min_m = med-4*std
    data = data[(data[:,1] < max_m) & (data[:,1] > min_m)]

    return data
