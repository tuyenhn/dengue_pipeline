#!/usr/bin/env bash

# Halt upon error
set -eoux pipefail
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Provide a name for the container here
NAME="grapevne-build"
TARGET_RULES="plotexportsandimports_target"

# Construct image
docker build \
    --platform linux/amd64 \
    --no-cache \
    --build-arg HOST_UID=$(id -u) \
    --build-arg TARGET_RULES="${TARGET_RULES}" \
    -t "${NAME}" \
    "${SCRIPT_DIR}"
