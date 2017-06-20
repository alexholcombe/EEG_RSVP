# EEG_RSVP

## Reading in raw data

The Psychopy program saves the data as a .txt file the name of which is the "subject name" and a date and time stamp. It is OK to use the participant's real name, because these files live in the rawData/ folder here which is never made public because it is listed in the .gitignore file.

The dataPreprocess/ directory contains the script that anonymizes the subject names using a cipher. Note that the key for the cipher is in a file in that directory on the local machine, but that file never becomes public because it is listed in the .gitignore file.

The dataPreprocess/loadAnonymizeSaveData.R file does what it says on the tin, putting the anonymised data in the dataAnonymized/ folder. The MOTcircular repo contains another version of this file that matches up datafiles with eyelink eyetracker files (no eyetracker done in this experiment, as yet).

## Issues

* Why are there NaN values for dat$responsePosRelative1?
* Should save the whole letter sequence, nost just the letters