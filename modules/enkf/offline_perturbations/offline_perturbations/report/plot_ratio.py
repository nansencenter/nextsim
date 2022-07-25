import scipy.io
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 

# import .mat data file
adjusted_ratio = scipy.io.loadmat('adjusted_ratio.mat')
# select keyword Ratio for the Ratio variable
Ratio = adjusted_ratio['Ratio']
# skip values when original wind files are missing.
Ratio = np.where(Ratio == 0, np.nan, Ratio)
# amplifer that U10_perturb = U10 + noise_U*amplifer
r = np.arange(0.5, 5.5, 0.5)
# set dates for plotting
start_time = dt.datetime(2019, 9, 4, 0, 0, 0)
dates = [(start_time + dt.timedelta(days = i*0.25)) for i in range(len(Ratio))]

# set font to time new roman
plt.rcParams["font.family"] = "serif"
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
# # set colormap for lines
cm = plt.get_cmap('tab10')

plt.figure(1, figsize=(16, 5))
plt.rcParams.update({'font.size': 16})

for i, rr in enumerate(r):
    plt.plot(dates, Ratio[:, i], label=f"r = {rr}", color=cm( i/len(r) ) )

plt.ylabel(r'ratio ($r$)', fontsize=24)
# x axis limit
plt.xlim([dt.date(2019, 6, 30), dt.date(2020, 10, 1)])
# locator etc.
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%d/%m/%Y"))
# location of legend
plt.legend(loc='center left')
# figure output
plt.savefig('ratio_curves.pdf')


# # %  --- plot temporal average ratio against r
plt.figure(2, figsize=(10, 8))
plt.rcParams.update({'font.size': 16})
plt.plot(r, np.nanmean(Ratio, 0), 'o-', markerfacecolor='None')
plt.xlabel(r'noise amplifer $r$', fontsize=24)
plt.ylabel(r'1-year averaged ratio $r$', fontsize=24)
plt.savefig("1-year_averaged_ratio.pdf")

