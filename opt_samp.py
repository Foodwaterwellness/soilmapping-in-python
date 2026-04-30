import numpy as np # matrix calculation library
import matplotlib.pyplot as plt # plotting library
from tqdm.notebook import trange, tqdm # progress bar library
import joblib as job # for multithreading on several cores at once

from scipy import stats, optimize, spatial # specialized math libraries
from sklearn import preprocessing, cluster # specialized machine learning libraries

import clhs # latin hypercube library

# Packaged algorithm for multithreading for cLHS
def _clhs_processing(covs, iterations, sample_size):
    result = clhs.clhs(covs, num_samples=sample_size, max_iterations=iterations, progress=True)
    return result["sample_indices"]

# Packaged algorithm for multithreading for FSCS
def _fscs_processing(covs, iterations, sample_size):
    covs_scaled = preprocessing.StandardScaler().fit_transform(covs)
    results = cluster.KMeans(n_clusters=sample_size, n_init=20, max_iter=iterations).fit(covs_scaled)

    distance = spatial.distance.cdist(results.cluster_centers_, covs_scaled)
    minimized = np.argmin(distance, axis=1)
    return minimized

# Use conditioned latin hypercube sampling for the analysis
def clhs_setup(cpus, clhs_iter, cseq, covs):
    parallel = job.Parallel(n_jobs=cpus, return_as="list", verbose=1)
    results = parallel(job.delayed(_clhs_processing)(covs, clhs_iter, i) for i in cseq)
    
    # Convert list output into numpy array
    clhs_results = []
    for array in results:
        clhs_results.append(array)
    return clhs_results

# Use feature space coverage sampling for the analysis
def fscs_setup(cpus, fscs_iter, cseq, covs):
    parallel = job.Parallel(n_jobs=cpus, return_as="list", verbose=1)
    results = parallel(job.delayed(_fscs_processing)(covs, fscs_iter, i) for i in cseq)
    
    # Convert list output into numpy array
    fscs_results = []
    for array in results:
        fscs_results.append(array)
    return fscs_results

# This standalone function calculates an estimate of the number of bins for the divergence calculations
def calculate_bins(covs):
    bins_FD = []
    bins_Scott = []
    for i in tqdm(np.arange(0, covs.shape[1])):
        q25, q75 = np.percentile(covs[:,i], [25,75])
        iqr = q75 - q25
        stdev = np.std(covs[:,i])
        
        bin_width_FD = 2*iqr/np.cbrt(covs.shape[0])
        bin_num_FD = int((np.max(covs[:,i]) - np.min(covs[:,i]))/bin_width_FD)
    
        bin_width_Scott = 3.49*stdev/np.cbrt(covs.shape[0])
        bin_num_Scott = int((np.max(covs[:,i]) - np.min(covs[:,i]))/bin_width_Scott)
        
        bins_FD.append(bin_num_FD)
        bins_Scott.append(bin_num_Scott)
    
    print(f"FD median bins: {np.median(bins_FD)}")
    print(f"Scott median bins: {np.median(bins_Scott)}")


    # Using the median value when calculating bins as stated in the paper
    bins = int(np.median(bins_FD))
    # bins = int(np.median(bins_Scott)) # Can return Scott calculation instead if desired.
    return bins

# Calculate the quantiles matrix from the covariates
def calculate_quantiles(bins, covs):
    quantiles = np.empty((bins+1, covs.shape[1]))
    for i in tqdm(np.arange(0, covs.shape[1])):
        quantiles[:,i] = np.linspace(np.min(covs[:,i]), np.max(covs[:,i]), bins+1)
    print(quantiles.shape)
    return quantiles

# Calculate the covariate matrix from the covariates
def calculate_covariates(bins, quantiles, covs):
    cov_matrix = np.ones((bins, covs.shape[1]))
    for i in tqdm(np.arange(0, covs.shape[0])):
        for j in np.arange(0, covs.shape[1]):
            dd = covs[i,j]
            for k in np.arange(0, bins):
                kl = quantiles[k,j]
                ku = quantiles[k+1,j]
                if dd >= kl and dd <= ku:
                    cov_matrix[k,j] = cov_matrix[k,j]+1
    print(cov_matrix.shape)
    return cov_matrix

# Main loop can use either cLHS or FSCS results
# Calculate the divergence metrics for the covariates
def calculate_divergences(setup_results, bins, covs, quantiles, cov_matrix):
    mean_storage = {}
    score_storage = {}
    for i in tqdm(setup_results):
        selected_samples = covs[i,:]
    
        h_matrix = np.ones((bins, covs.shape[1]))
        for i in np.arange(0, selected_samples.shape[0]):
            for j in np.arange(0, selected_samples.shape[1]):
                dd = selected_samples[i,j]
                for k in np.arange(0, bins):
                    kl = quantiles[k,j]
                    ku = quantiles[k+1,j]
                    if dd >= kl and dd <= ku:
                        h_matrix[k,j] = h_matrix[k,j]+1
    
        kl_storage = np.empty((covs.shape[1], 1))
        for i in np.arange(0, covs.shape[1]):
            kl_storage[i] = np.sum(stats.entropy(np.divide(cov_matrix[:,i], np.sum(cov_matrix[:,i])), np.divide(h_matrix[:,i], np.sum(h_matrix[:,i])), base=2))
        kldiv_result = np.mean(kl_storage)
    
        js_storage = np.empty((covs.shape[1], 1))
        for i in np.arange(0, covs.shape[1]):
            left_dist = np.divide(cov_matrix[:,i], np.sum(cov_matrix[:,i]))
            right_dist = np.divide(h_matrix[:,i], np.sum(h_matrix[:,i]))
            middle_dist = np.multiply(0.5, np.add(left_dist, right_dist))
            js_storage[i] = 0.5*np.sum(stats.entropy(left_dist, middle_dist, base=2)) + 0.5*np.sum(stats.entropy(right_dist, middle_dist, base=2))
        jsdiv_result = np.mean(js_storage)
        jsdist_result = np.mean(np.sqrt(js_storage))
    
        if selected_samples.shape[0] in mean_storage:
            mean_storage[selected_samples.shape[0]].append([kldiv_result, jsdiv_result, jsdist_result])
            score_storage[selected_samples.shape[0]].append([kl_storage, js_storage, np.sqrt(js_storage)])
        else:
            mean_storage[selected_samples.shape[0]] = [[kldiv_result, jsdiv_result, jsdist_result]]
            score_storage[selected_samples.shape[0]] = [[kl_storage, js_storage, np.sqrt(js_storage)]]

    # Repackage and summarize the divergence values into lists for easier curve fitting and plotting
    kldiv_storage = []
    jsdiv_storage = []
    jsdist_storage = []
    
    for c in mean_storage:
        kldiv = []
        jsdiv = []
        jsdist = []
        for array in mean_storage[c]:
            kldiv.append(array[0])
            jsdiv.append(array[1])
            jsdist.append(array[2])
        kldiv_storage.append(np.mean(kldiv))
        jsdiv_storage.append(np.mean(jsdiv))
        jsdist_storage.append(np.mean(jsdist))

    return kldiv_storage, jsdiv_storage, jsdist_storage

# Model function to fit the exponential decay of the divergence
def _exp_func(x, a, b, c):
    return a + b*np.exp(-c*x)

# Fitting the divergence function for each of the metrics and then plotting the results
def plot_divergences(cseq, conf, kldiv_storage, jsdiv_storage, jsdist_storage):
    nseq = []
    for c in cseq:
        if c not in nseq:
            nseq.append(c)
    
    fit_kldiv, fit_cov = optimize.curve_fit(_exp_func, nseq, kldiv_storage, p0=[0.01, np.max(kldiv_storage)*0.9, 0.01])
    fit_jsdiv, fit_cov = optimize.curve_fit(_exp_func, nseq, jsdiv_storage, p0=[0.01, np.max(jsdiv_storage)*0.9, 0.01])
    fit_jsdist, fit_cov = optimize.curve_fit(_exp_func, nseq, jsdist_storage, p0=[0.01, np.max(jsdist_storage)*0.9, 0.01])
    
    test_x = np.arange(0, np.max(cseq))
    
    y_kldiv = fit_kldiv[0] + fit_kldiv[1]*np.exp(-fit_kldiv[2]*test_x)
    y_jsdiv = fit_jsdiv[0] + fit_jsdiv[1]*np.exp(-fit_jsdiv[2]*test_x)
    y_jsdist = fit_jsdist[0] + fit_jsdist[1]*np.exp(-fit_jsdist[2]*test_x)

    plt.figure()
    plt.plot(nseq, kldiv_storage, "r+")
    plt.plot(nseq, jsdiv_storage, "r+")
    plt.plot(nseq, jsdist_storage, "r+")
    plt.plot(test_x, y_kldiv, label="KL Div")
    plt.plot(test_x, y_jsdiv, label="JS Div")
    plt.plot(test_x, y_jsdist, label="JS Dist")
    plt.legend()
    plt.ylabel("Metric value")
    plt.xlabel("Number of samples")
    plt.title("Exponential Decay of Three Metrics")

    # Calculate a simple CDF for the exponential decay and find the value where the confidence value is reached
    norm_kldiv = 1 - (y_kldiv - np.min(y_kldiv))/(np.max(y_kldiv) - np.min(y_kldiv))
    norm_jsdiv = 1 - (y_jsdiv - np.min(y_jsdiv))/(np.max(y_jsdiv) - np.min(y_jsdiv))
    norm_jsdist = 1 - (y_jsdist - np.min(y_jsdist))/(np.max(y_jsdist) - np.min(y_jsdist))

    plt.figure()
    plt.plot(test_x, norm_kldiv, label="KL Div")
    plt.plot(np.where(norm_kldiv > conf)[0][0], norm_kldiv[np.where(norm_kldiv > conf)[0][0]], "+r")
    plt.plot(test_x, norm_jsdiv, label="JS Div")
    plt.plot(np.where(norm_jsdiv > conf)[0][0], norm_jsdiv[np.where(norm_jsdiv > conf)[0][0]], "+r")
    plt.plot(test_x, norm_jsdist, label="JS Dist")
    plt.plot(np.where(norm_jsdist > conf)[0][0], norm_jsdist[np.where(norm_jsdist > conf)[0][0]], "+r")
    plt.grid()
    plt.legend()
    plt.ylabel("CDF value")
    plt.xlabel("Number of samples")
    plt.title(f"CDF Functions with {conf} Cutoff")
    print(f"KL div optimal n: {np.where(norm_kldiv > conf)[0][0]}")
    print(f"JS div optimal n: {np.where(norm_jsdiv > conf)[0][0]}")
    print(f"JS dist optimal n: {np.where(norm_jsdist > conf)[0][0]}")