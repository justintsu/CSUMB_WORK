import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

from IPython import embed as shell
import sys

# Define a function to calculate cumulative stats
def accum_data(df,typ='Mean_flux'):
    # Requires a dataframe containing N2O data in ascending order
    #  date from one site and one treatment
    # Returns a value for cumulative N2O based on equation:
    #  C = 0.5 * (D_a + D_b) * (B - A)
    # (B-A)
    time_diff = [-td/np.timedelta64(1,'D') for td in df.index.diff(periods=-1).dropna()]

    # (D_a + D_b)
    flux_sum = df[typ].rolling(window=2).sum().dropna().values

    return [0.5*f*t for f,t in zip(flux_sum, time_diff)]

def main():
    parser = argparse.ArgumentParser(description="A plotting script that takes in a CSV file with site, treatment and flux data.")

    parser.add_argument("-f", type=str, required=True, help="Pathway to CSV Data File")
    parser.add_argument("-s", type=str, nargs="+", required=False, help="User Specified Sites to Plot (All by Default)")
    parser.add_argument("-t", type=str, nargs="+", required=False, help="User Specified Treatments to Plot (All by Default)")

    args = parser.parse_args()
    # Specify CSV format input data file
    f_input = args.f

    # Import CSV format data using pandas read_csv function
    df = pd.read_csv(f_input)

    shell()
    # Drop the garbage unnamed column
    df.drop([df.columns[0]], axis=1, inplace=True)

    # Convert the date column to a datetime object
    df['Date'] = pd.to_datetime(df['Date'])

    # Reset the index of the dataframe to date
    df.set_index(df.columns[0], inplace=True)

    # Isolate only the nitrous oxide species data
    df = df[df['Species'] == 'N2O']

    # Initialize figure and axes
    fig = plt.figure(figsize=(12,6))
    fig.suptitle('Cumulative Flux and Error by Treatment and Site')

    site_labels, flux_values, err_values = {}, {}, {}

    # Determine which treatments to loop through based on user specs
    # If '-t' was not supplied, use all
    shell()
    sys.exit()
    if args.t:
        treatments = args.t
    else:
        treatments = df.Treatment.unique()

    # Loop through Treatment Types
    for treatment in treatments:
        #Check to see if treatment is in the dataframe, if not then skip loop iteration
        if treatment not in df['Treatment'].unique(): continue

        # Plotting stuff
        site_labels[treatment] = []
        flux_values[treatment] = []
        err_values[treatment] = []

        # Determine which Sites to loop through based on user specs
        # If '-s' was not supplied, use all
        if args.s:
            sites = args.s
        else:
            sites = df[df['Treatment'] == treatment].Site.unique()

        # Loop through Sites within that Treatment Type
        for site in sites:
            # Check to see if site is in the dataframe, if not then skip loop iteration
            if site not in df['Site'].unique(): continue
            data = df[(df['Treatment'] == treatment) & (df['Site'] == site)].sort_index(ascending=True)
            cumulative_fluxes = accum_data(data, 'Mean_flux')
            # Convert fluxes from mg/m2 to kg/hectare
            cumulative_fluxes = [cumulative_flux *  (1.e4 / 1.e6) for cumulative_flux in cumulative_fluxes]

            # Add up all fluxes to get total flux
            cumulative_flux = np.sum(cumulative_fluxes)

            # Compute standard deviation of cumulative fluxes
            cumulative_flux_se = np.std(cumulative_fluxes) / np.sqrt(len(cumulative_fluxes))
            #
            
            start = data.index.min().strftime('%m/%d/%Y')
            end = data.index.max().strftime('%m/%d/%Y')
            print('Site {4}, Treatment {5}:\n' \
                  '{2} - {3}\n' \
                  '_______________________________ \n' \
                  'Cumulative Flux: {0} \n' \
                  'Standard Error: {1}\n\n\n'.format(cumulative_flux, cumulative_flux_se, start, end, site, treatment))

            # Plotting stuff
            site_labels[treatment].append(f'{site}\n{start} - {end}')
            flux_values[treatment].append(cumulative_flux)
            err_values[treatment].append(cumulative_flux_se)

    # Bar positions
    num_treatments = len(site_labels.keys())
    treatment_indices = np.arange(num_treatments)
    width = 0.2 * num_treatments

    # Define a finite list of colors (number of treatments cannot exceed length of this list)
    color_list = ['blue', 'red', 'cyan', 'yellow', 'magenta']



    #Combine the site_label data with flux_values and err_values
    flux = {}
    se = {}
    for treatment in site_labels.keys():
        flux[treatment] = {}
        se[treatment] = {}
        for i,site in enumerate(site_labels[treatment]):
            flux[treatment][site] = flux_values[treatment][i]
            se[treatment][site] = err_values[treatment][i]

    # Make a master site list that includes all sites from all treatments
    site_list = np.unique(sum(list(site_labels.values()), []))
    # Make a width array based on the treatments
    rel_bar_pos = np.linspace(-1. * width / 2., width / 2., len(flux.keys()))

    for treatment_number, treatment in enumerate(flux.keys()):
        # treatment_number will start at 0 and go to N-1 # of treatments
        #  this value can be used as a proxy for the index of treatment color
        #  and as a proxy for the bar positioning in each site cluster
        bar_pos = rel_bar_pos[treatment_number]
        for site_tick,site in enumerate(site_list):
            # Check if the site in the master site_list is in your flux key
            if site not in flux[treatment].keys(): continue
            plt.bar(site_tick+bar_pos, flux[treatment][site], width, label=f'{treatment}',
                    color=color_list[treatment_number], alpha=0.7, edgecolor='black', linewidth=1.5)
            plt.errorbar(site_tick+bar_pos, flux[treatment][site], yerr=se[treatment][site],
                         fmt="o", color="black")

    # Formatting
    plt.xticks(np.arange(len(site_list)),site_list,rotation=45.)
    plt.ylabel('N2O kg/ha')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.tight_layout()
    plt.savefig('flux_and_err.png')
    plt.show()


if __name__ == "__main__":
    main()
