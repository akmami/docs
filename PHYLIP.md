# PHYLIP

This documentation provides information related to the PHYLIP tool. The latest version of PHYLIP is 3.697.

**Note:**  
- This guide is tailored for Linux and Unix systems.  
- Mac users should refer to relevant installation guides for their platform.  
- Windows users may cry while writing this in their diary.

## Installation

To install PHYLIP, follow these steps:

1. Download the source code from the following link:  
   [PHYLIP 3.697 Source Code](https://phylipweb.github.io/phylip/download/phylip-3.697.tar.gz)

2. Extract the downloaded files by running:  
   ```sh
   tar -zxvf phylip-3.697.tar.gz
   cd phylip-3.697/src
   ```

3. Create a `Makefile` by copying the existing template:  
   ```sh
   cat Makefile.unx > Makefile
   ```

4. Modify **line 94** in the `Makefile` to resolve potential compilation issues. Replace the `CFLAGS` definition with the following:  
   ```txt
   CFLAGS = -zmuldefs
   ```

5. Compile the source code:  
   ```sh
   make install
   ```

## Running PHYLIP

The compiled program executables are located in the `exe` directory. Navigate to this directory to run the desired tool.

**Important Notes:**  
- PHYLIP does not include a graphical user interface (GUI).  
- All operations must be performed through the terminal.

For additional information regarding the available commands, refer to the official [PHYLIP documentation](https://phylipweb.github.io/phylip/phylip.html).
