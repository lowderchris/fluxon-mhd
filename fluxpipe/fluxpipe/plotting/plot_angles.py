import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('qt5Agg')
# Specify the path to your CSV file
csv_file_path = 'seq_ass_theta.csv'

# Specify your column names
column_names = ['Fid', 'ass', 'theta', 'speed']

# Read the CSV file
df = pd.read_csv(csv_file_path, header=None, names=column_names, sep=", ")

# Display the DataFrame
# print(df)


# Plotting
df.plot(kind='scatter', x='theta', y='ass', c='speed', alpha=0.5, edgecolor="none")
plt.title(r'Scatter Plot of $A_{ss}$ vs $\theta_b$, ' + str(len(df)) + " points" )
plt.xlabel('Edge Angle')
plt.xlim((0, 80))
plt.ylim((0, 100))
plt.ylabel('Expansion Factor')
plt.axvline(6, ls="--")
plt.axvline(10)
# plt.show()

# Histogram
plt.figure(figsize=(10, 6))  # Create a new figure for the histogram
plt.hist(df['speed'], bins=10, alpha=0.7, edgecolor='black')
plt.title('Histogram of Speeds')
plt.xlabel('Speed')
plt.ylabel('Frequency')
plt.show()