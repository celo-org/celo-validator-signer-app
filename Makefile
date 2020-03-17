include $(BOLOS_SDK)/Makefile.defines

APPNAME = "Celo Validator"
APPVERSION = 1.0.0
APP_LOAD_PARAMS = --appFlags 0x200 $(COMMON_LOAD_PARAMS)

all: ledger-app

ledger-app: dockerimage bls-embedded
	docker run -v ${PWD}:/code celo-org/validator-signer "cd ledger-app && make all"

fp-opt: 
	docker run -v ${PWD}:/code celo-org/validator-signer "./compiling.sh && cd bls-embedded/bls12_377 && make"

bls-embedded: fp-opt
	docker run -v ${PWD}:/code celo-org/validator-signer "cd bls-embedded/bls && RUSTFLAGS=\"-L /code/bls-embedded/bls12_377\" cargo +nightly build --release --target thumbv6m-none-eabi"

clean:
	docker run -v ${PWD}:/code celo-org/validator-signer "cd ledger-app && make clean && cd ../bls-embedded/bls && cargo clean"

dockerimage:
	docker build -t 'celo-org/validator-signer' .

load: all
	TARGET_ID=0x33000004
	cd ledger-app && python -m ledgerblue.loadApp $(APP_LOAD_PARAMS)

delete:
	cd ledger-app && python -m ledgerblue.deleteApp $(COMMON_DELETE_PARAMS)
