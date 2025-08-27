# Gungnir
Guess till correct: Gungnir codec enabling low redundancy DNA storage through substantial computing power

## Introduction
Gungnir is a DNA storage codec system that supports adjustable configurations to maintain a balance between information density and error tolerance. Powered by the computing power, Gungnir catches up with the performance of traditional methods with the use of only half the DNA bases. 

This repository provides a complete toolkit for running sample test, including encoding/decoding algorithms, simulation tools, and practical examples. 

## Content
- [Introduction](#introduction)
- [Installation](#installation)
  - [Install go 1.23](#1install-go-123)
  - [Get Gungnir source code](#2get-gungnir-source-code)
- [Quick Start](#quick-start)
  - [DNA Encoding](#1-dna-encoding)
  - [Add Noise](#2-add-noise)
  - [DNA Decoding](3-dna-decoding)
  - [Reconstruction](#4-reconstruction)
- [Customized Usage](#customized-usage)

## Files Tree Diagram
```
Gungnir
├── error_pattern
│   └── ...                                  # Error handling patterns
├── examples
│   └── main.go                              # Example usage code
├── files
│   └── The Ugly Duckling                    # Sample test file
├── tools
│   ├── decode.go                            # DNA decoding
│   ├── decode_three.go                      # Ternary DNA decoding
│   ├── distance.go                          # Distance calculation
│   ├── encode.go                            # DNA encoding
│   ├── exclude.go                           # Invalid motifs
│   ├── hash.go                              # Hash functions
│   ├── params.go                            # Parameters
│   ├── readfile.go                          # File reading functions
│   └── simulation.go                        # Error simulation
├── .gitignore                               # Git ignore
├── LICENSE                                  # Project license
├── README.md                                # Description file
├── go.mod                                   # Go module
└── go.sum                                   # Dependency checksums
```
## Installation
### 1.Install go 1.23
Download and install Golang from the official website (https://golang.org/dl/).
```
# Download the latest Golang version 1.23 by visiting the official website (https://golang.org/dl/) and, 
# copying the download link for the Linux tarball.
# An example is shown below:
wget https://golang.org/dl/go1.23.10.linux-amd64.tar.gz

# Extract the downloaded tarball to your preferred local directory. In this example, we'll use `$HOME/.local`:
mkdir -p $HOME/.local
tar -xvzf go1.23.10.linux-amd64.tar.gz -C $HOME/.local

# Remove the tarball after extraction
rm go1.23.10.linux-amd64.tar.gz

# Set up your Go workspace and environment variables
## Create the required directory structure:
mkdir -p $HOME/go/{bin,src,pkg}

## add link to bashrc or .profile
## add the GOPATH,GOROOT to your `~/.bashrc` or `~/.profile`
echo 'export GOPATH=$HOME/go
export GOROOT=$HOME/.local/go
export PATH=$PATH:$GOROOT/bin:$GOPATH/bin' >> ~/.bashrc
source ~/.bashrc

## Verify the installation
## You will get "go version go1.23.10 linux/amd64" if installed successfully
go version
```
### 2.Get Gungnir source code
You can clone this repo as following:
```
mkdir Gungnir_RootFolder
cd Gungnir_RootFolder
git clone git@github.com:HKU-BAL/Gungnir.git
cd Gungnir

# $Gungnir_DIR is path of Gungnir
Gungnir_DIR=$(pwd)
```
## Quick Start
This section demonstrates the basic workflow of our DNA storage system. The following commands will encode a file into DNA sequences, simulate sequencing errors, decode the noisy sequences, and reconstruct the original file.
### 1. DNA Encoding
Convert your digital file into DNA sequences (A, T, C, G). This process will generate an *Origin* file inside the *Outcome* directory.
```
go run main.go -action Encode -input "../files/Summer Flowers" -output "../Outcome"
```
Default parameters apply the Gungnir method at 0.8 density, you can modify these by applying optional parameters.
### 2. Add Noise
Simulate DNA sequencing errors aiming at testing the robustness of the codec.
```
go run main.go -action AddNoise -output "../Outcome"
```
This adds 1% substitution, 1% insertion, and 1% deletion errors by default.

### 3. DNA Decoding
Decode the noisy DNA sequences back to digital data with error correction.
```
go run main.go -action Decode -output "../Outcome"
```
The decoder automatically corrects errors using adaptive edit distance algorithms.
### 4. Reconstruction
Rebuild the original file from the decoded data.
```
go run main.go -action Reconstruction -output "../Outcome"
```
After this step, your original file will be recovered in *output* file inside the *Outcome* directory.

## Customized Usage
Our codec supports various modifiable parameters to optimize specific use case:

