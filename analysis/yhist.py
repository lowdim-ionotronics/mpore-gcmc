# %%
import numpy as np
import pandas as pd
import warnings 
warnings.filterwarnings("ignore") 
import sys

# %%
def read_mctrj(fn):

    data_df = pd.DataFrame(columns=['q', 'x', 'y', 'z'])

    data = pd.read_csv(fn, delimiter='\s+', skiprows=0, on_bad_lines='skip', names=['q', 'x', 'y', 'z'])

    print(data.head(10))
    cut1 = np.array(data['x'].isnull())
    cut1id = np.array(np.where(cut1==True)[0])

    traj = data[['q','x','y','z']].to_numpy(copy=True)
    frame_traj_values = []

    for i in np.arange(0, cut1id.size-1):
        iframe_traj_values = traj[cut1id[i]+1:cut1id[i+1], :].astype(float)
        iframe_traj_values = iframe_traj_values.astype(float)
        iframe_traj_values = iframe_traj_values.reshape(-1,4)
        print(iframe_traj_values.shape[0])
        frame_traj_values.append(iframe_traj_values) 
        
    return frame_traj_values


# %%
def write_yhist(values, frames, fn):

    z_bin_edges = np.arange(-10.0,10.0,0.05) 
    z_bin_centers = (z_bin_edges[1:] + z_bin_edges[:-1])/2.0
    z_hist_list = np.empty([0])

    for i in np.arange(len(values)-frames, len(values)):

        value = values[i]
        z_locs = value[:, 3]
        frame_hist, _ = np.histogram(z_locs, z_bin_edges)
        z_hist_list = np.append(z_hist_list, frame_hist)

    z_hist_list = z_hist_list.reshape(frames, -1)
    z_hist_mean = np.mean(z_hist_list, axis=0)

    np.savetxt(fn, np.c_[z_bin_centers, z_hist_mean])

# %%
inname = sys.argv[1]
outname = sys.argv[2]
frames = int(sys.argv[3])

# inname = "../asymmetric_equil_1.0r_1l_3.5M_0.5p.lammpstrj"
# outname = "test.q"

values = read_mctrj(inname)
write_yhist(values, frames, outname)



