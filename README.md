# SEC_C

Super-Efficient Cross-Correlation (SEC-C): a fast matched filtering code for seismic waveforms, optimized for use on desktop computers

SEC-C is an ongoing effort for accelerating the speed of cross-correlation analysis for seismological applications. We have a manuscript in the review. If you are using SEC-C for your work/research, please cite the following paper:

Nader Shakibay Senobari, Gareth J. Funning, Eamonn Keogh, Yan Zhu, Chin‐Chia Michael Yeh, Zachary Zimmerman, Abdullah Mueen; Super‐Efficient Cross‐Correlation (SEC‐C): A Fast Matched Filtering Code Suitable for Desktop Computers. Seismological Research Letters doi: https://doi.org/10.1785/0220180122

SEC-C is written in MATLAB, but a python version is also provided. The current python version is slower than MATLAB version, but the work is in progress.

SEC-C has two main branches: SEC-C for continues data (i.e., template matching or matched filtering) and SEC-C for individual waveforms (i.e., pairwise cross-correlation). Former is already uploaded, and we are going to upload the individual case soon.

A toy example of performing template matching that includes, retrieving, prepossessing, performing template matching using SEC-C and postprocessing results for Mt St Helens seismicity is now included (Mt_St_Helens_example.m).
