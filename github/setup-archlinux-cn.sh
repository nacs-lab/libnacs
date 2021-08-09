#!/bin/bash

cat >> /etc/pacman.conf <<EOF
[archlinuxcn]
Server = https://repo.archlinuxcn.org/\$arch
EOF

pacman -Sy --noconfirm archlinuxcn-keyring
pacman -Sy --noconfirm devtools-cn
