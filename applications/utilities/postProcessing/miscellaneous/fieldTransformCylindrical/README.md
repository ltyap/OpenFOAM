# README #

This is an adaptation of the exisitng OpenFOAM utility called `postChannelFlow` for OpenFOAM v9, which allows one to average the solution of a channel flow simulation over the stream- and spanwise directions. The original utility can be found [here](https://github.com/timofeymukha/postChannelFlow)
To run the utility just write postChannelFlow in the command-line.
The utility expectes a file called postChannelDict in the constant directory.
A sample dictionary can be found in the repository.

Some differences with the `postChannel`
- Averages all the fields you have in the time directory.
- Averages data on the wall patches.
- Does *not* average across the channel centerline, so you get the full profile across the channel.

**This offering is not approved or endorsed by OpenCFD Limited, producer
and distributor of the OpenFOAM software and owner of the OPENFOAM®  and
OpenCFD®  trade marks.**