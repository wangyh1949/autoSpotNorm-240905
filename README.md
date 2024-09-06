autoSpotNorm combines the result from oufti & uTrack.
it calculates the normalized position of spots and diffusion properties of tracks. 

More specifically, autoSpotNorm.m automatically run spotNorm.m for all cell meshes & tracking movies in the folder, and combine them into the same tracksFinal variable. It will then calcuate quantities for analysis (e.g. MSD, Diff, alpha). Eventually, these quantities are converted into matrices (arrays) which are convenient for future plotting
