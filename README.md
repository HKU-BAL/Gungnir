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
- [Data Preprocessing](#data-preprocessing)
- [Quick Start](#quick-start)
  - [DNA Encoding](#dna-encoding)
  - [DNA Decoding](#dna-decoding)
  - [Distance Calculation](#distance-calculation)
  - [Error Simulation](#error-simulation)

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






