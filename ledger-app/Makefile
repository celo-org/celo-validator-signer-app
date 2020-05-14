#*******************************************************************************
#   Ledger Blue
#   (c) 2016 Ledger
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#*******************************************************************************


ifeq ($(BOLOS_SDK),)
$(error BOLOS_SDK is not set)
endif

include $(BOLOS_SDK)/Makefile.defines

SIGNER_LIB := ../bls-embedded/bls/target/thumbv6m-none-eabi/release/libbls_embedded.a ../bls-embedded/bls12_377/libfpc.a
SIGNER_LIB_INCLUDES := ../bls-embedded/bls ../bls-embedded/bls12_377/src

# Main app configuration

APPNAME = "Celo Signer"
APPVERSION = 1.0.0
APP_LOAD_PARAMS = --appFlags 0x200 $(COMMON_LOAD_PARAMS)

# Build configuration

APP_SOURCE_PATH += src
SDK_SOURCE_PATH += lib_stusb lib_stusb_impl

DEFINES += APPVERSION=\"$(APPVERSION)\"

DEFINES += OS_IO_SEPROXYHAL IO_SEPROXYHAL_BUFFER_SIZE_B=128
DEFINES += HAVE_BAGL HAVE_SPRINTF
DEFINES += PRINTF\(...\)=

DEFINES += HAVE_IO_USB HAVE_L4_USBLIB IO_USB_MAX_ENDPOINTS=7 IO_HID_EP_LENGTH=64 HAVE_USB_APDU

DEFINES += CX_COMPLIANCE_141

# Specific to nano X
DEFINES       += HAVE_BLE BLE_COMMAND_TIMEOUT_MS=2000
DEFINES       += HAVE_BLE_APDU # basic ledger apdu transport over BLE

# Specific to nano X
DEFINES       += HAVE_GLO096
DEFINES       += HAVE_BAGL BAGL_WIDTH=128 BAGL_HEIGHT=64
DEFINES       += HAVE_BAGL_ELLIPSIS # long label truncation feature
DEFINES       += HAVE_BAGL_FONT_OPEN_SANS_REGULAR_11PX
DEFINES       += HAVE_BAGL_FONT_OPEN_SANS_EXTRABOLD_11PX
DEFINES       += HAVE_BAGL_FONT_OPEN_SANS_LIGHT_16PX

# Compiler, assembler, and linker

ifneq ($(BOLOS_ENV),)
$(info BOLOS_ENV=$(BOLOS_ENV))
CLANGPATH := $(BOLOS_ENV)/clang-arm-fropi/bin/
GCCPATH := $(BOLOS_ENV)/gcc-arm-none-eabi/bin/
else
$(info BOLOS_ENV is not set: falling back to CLANGPATH and GCCPATH)
endif

ifeq ($(CLANGPATH),)
$(info CLANGPATH is not set: clang will be used from PATH)
endif

ifeq ($(GCCPATH),)
$(info GCCPATH is not set: arm-none-eabi-* will be used from PATH)
endif

CC := $(CLANGPATH)clang
CFLAGS += -O3 -Os

AS := $(GCCPATH)arm-none-eabi-gcc
AFLAGS +=

LD := $(GCCPATH)arm-none-eabi-gcc
LDFLAGS += -O3 -Os -L/opt/bolos-env/gcc-arm-none-eabi/arm-none-eabi/lib/thumb/
LDLIBS += -lm -lc

# Specific to nano X
SDK_SOURCE_PATH  += lib_blewbxx lib_blewbxx_impl
SDK_SOURCE_PATH  += lib_ux
DEFINES          += HAVE_UX_FLOW

include $(BOLOS_SDK)/Makefile.glyphs

# Main rules

all: default

# Import generic rules from the SDK

include $(BOLOS_SDK)/Makefile.rules

OBJECT_FILES += $(SIGNER_LIB)
INCLUDES_PATH += $(SIGNER_LIB_INCLUDES)
