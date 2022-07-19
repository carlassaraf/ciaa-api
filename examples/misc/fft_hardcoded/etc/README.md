# etc

This directory contains the [lpc4337.cfg](openocd/lpc4337.cfg) file to debug the CIAA with the MCUXpresso.

First, copy this directory into your new project's root directory.

Then, set up a new `GDB OpenOCD Debugging` configuration and within the debugging tab add the next line to config options:

```
-f etc/openocd/lpc4337.cfg
```

Now the debugging is available.