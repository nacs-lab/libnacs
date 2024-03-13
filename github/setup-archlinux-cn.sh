#!/bin/bash

cat >> /etc/pacman.conf <<EOF
[archlinuxcn]
Server = https://repo.archlinuxcn.org/\$arch
EOF

# https://github.com/archlinuxcn/repo/issues/3557
# https://www.archlinuxcn.org/archlinuxcn-keyring-manually-trust-farseerfc-key/
pacman-key --lsign-key "farseerfc@archlinux.org"

pacman -Sy --noconfirm archlinuxcn-keyring
pacman -Sy --noconfirm devtools-cn
