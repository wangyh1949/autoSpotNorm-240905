## autoSpotNorm

autoSpotNorm calculates the normalized position of spots and diffusion properties of tracks. It combines the output from oufti & uTrack and carry out further analysis.

More specifically, autoSpotNorm.m automatically run spotNorm.m for all cell meshes & tracking movies in the folder, and combine them into the same tracksFinal variable. It will then calculate quantities for analysis (e.g. MSD, Diff, alpha). Eventually, these quantities are converted into matrices (arrays) which are convenient for future plotting.

For calculating the normalized position of spots within cells, check out spotNorm_yh.m for more details.

combineTF.m combines results from multiple days of SPT experiments