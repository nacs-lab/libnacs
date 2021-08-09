#!/bin/bash

# Adapted from FFY00/build-arch-package

# Initialize keyring
pacman-key --init
pacman-key --populate archlinux

pacman -Syu --noconfirm devtools
systemd-machine-id-setup # devtools requires a machine ID
sed -i "s|MAKEFLAGS=.*|MAKEFLAGS=-j$(nproc)|" /etc/makepkg.conf
useradd -m user # makepkg needs to be run with non-root
