import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    parser.add_argument("-o", action="store_true", required=False, default=False, help="Boolean to control cumulative flux CSV output (False by Default)")

    args = parser.parse_args()
    if args.o:
        print_csv = args.o
    else:
        print_csv = False

    # Specify CSV format input data file
    f_input = args.f

    # Import CSV format data using pandas read_csv function
    df = pd.read_csv(f_input)

    # Check for the existence of Date, Site, Treatment, Species,
    #   Reporting_units, and Mean_flux Columns
    column_check = ['Date', 'Site', 'Treatment', 'Species', 'Reporting_units', 'Mean_flux']
    for col in column_check:
        if col not in df.columns:
            print(f'{col} column does not exist in {f_input}. Exiting')
            sys.exit()

    # Drop the garbage unnamed column
    df.drop([df.columns[0]], axis=1, inplace=True)

    # Convert the date column to a datetime object
    df['Date'] = pd.to_datetime(df['Date'])

    # Reset the index of the dataframe to date
    df.set_index('Date', inplace=True)

    # Isolate only the nitrous oxide species data
    df = df[df['Species'] == 'N2O']


    site_labels, flux_values, err_values = {}, {}, {}

    # Determine which treatments to loop through based on user specs
    # If '-t' was not supplied, use all
    if args.t:
        treatments = args.t
    else:
        treatments = df.Treatment.unique()

    # Loop through Treatment Types
    printout_list = []
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
            cumulative_flux_se = np.std(cumulative_fluxes, ddof=1) / np.sqrt(len(cumulative_fluxes))
            #
        
            start = data.index.min().strftime('%m/%d/%Y')
            end = data.index.max().strftime('%m/%d/%Y')
            if not print_csv:
                #If we aren't writing a file, then print to screen
                print('Site {4}, Treatment {5}:\n' \
                      '{2} - {3}\n' \
                      '_______________________________ \n' \
                      'Cumulative Flux: {0} \n' \
                      'Standard Error: {1}\n\n\n'.format(cumulative_flux, cumulative_flux_se, start, end, site, treatment))
            else:
                printout_list.append([treatment, site, cumulative_flux, cumulative_flux_se])
            # Plotting stuff
            site_labels[treatment].append(f'{site}\n{start} - {end}')
            flux_values[treatment].append(cumulative_flux)
            err_values[treatment].append(cumulative_flux_se)

    #Create csv dataframe and export
    if print_csv:
        df_print = pd.DataFrame(printout_list)
        start = df.index.min().strftime('%Y%m%d')
        end = df.index.max().strftime('%Y%m%d')
        df_print.rename(columns={0:"Treatment", 1:"Site", 2:"Cumulative Flux", 3:"SE"}, inplace=True)
        df_print.to_csv(f'flux_and_err_{start}_{end}.csv', index=False)
    # Reorganize data
    num_treatments = len(site_labels.keys())
    treatment_indices = np.arange(num_treatments)
    # Make a master site list that includes all sites from all treatments
    site_list = np.unique(sum(list(site_labels.values()), []))
    #Combine the site_label data with flux_values and err_values
    flux = {}
    se = {}
    for treatment in site_labels.keys():
        flux[treatment] = {}
        se[treatment] = {}
        for i,site in enumerate(site_labels[treatment]):
            flux[treatment][site] = flux_values[treatment][i]
            se[treatment][site] = err_values[treatment][i]

    # Extract all uniquely named Sites
    unique_sites = sorted({site for treatment in flux.values() for site in treatment.keys()})

    # Extract all uniquely named Treatments
    unique_treatments = sorted(flux.keys())

    # Assign a unique color to each unique Treatment
    colors = sns.color_palette("tab10", len(unique_treatments))
    treatment_colors = dict(zip(unique_treatments, colors))

    # Prepare data for each bar
    site_index = {site: i for i,site in enumerate(unique_sites)}
    bar_width = 0.8/len(unique_treatments)
    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle('Cumulative Flux and Error by Treatment and Site')

    for i, (treatment, site_dict) in enumerate(flux.items()): #determine the Treatment bar we're working on
        for site, value in site_dict.items(): #determine which Site cluster this Treatment bar goes in
            # Adjust the position of the Treatment bar in the Site cluster
            x_pos = site_index[site] + (i - len(unique_treatments)/2) * bar_width
            error = se[treatment][site] if treatment in se and site in se[treatment] else 0
            ax.bar(x_pos, value, width=bar_width, color=treatment_colors[treatment], 
                   yerr=error, capsize=5)
               
    for treatment, color in treatment_colors.items():
        ax.bar(0, 0, color=color, label=treatment)

    ax.set_xticks(range(len(unique_sites)))
    ax.set_xticklabels(unique_sites, rotation=45, ha="right")
    ax.set_ylabel("N2O kg/ha")
    ax.legend(title="Treatments", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    plt.savefig('flux_and_err.png')
    plt.show()


if __name__ == "__main__":
    main()
