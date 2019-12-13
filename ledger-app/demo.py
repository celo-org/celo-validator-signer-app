#!/usr/bin/env python
#*******************************************************************************
#*   Ledger Blue
#*   (c) 2016 Ledger
#*
#*  Licensed under the Apache License, Version 2.0 (the "License");
#*  you may not use this file except in compliance with the License.
#*  You may obtain a copy of the License at
#*
#*      http://www.apache.org/licenses/LICENSE-2.0
#*
#*  Unless required by applicable law or agreed to in writing, software
#*  distributed under the License is distributed on an "AS IS" BASIS,
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#*  See the License for the specific language governing permissions and
#*  limitations under the License.
#********************************************************************************
from ledgerblue.comm import getDongle
from ledgerblue.commException import CommException
from secp256k1 import PublicKey
import time
import array

textToSign = str(bytearray([0, 124, 213, 66, 221, 251, 203, 196, 218, 98, 62, 47, 173, 180, 24, 40, 5, 115, 10, 177, 169, 95, 211, 218, 196, 160, 222, 2, 93, 129, 247, 93, 212, 138, 181, 114, 101, 206, 8, 244, 209, 107, 70, 47, 8, 231, 2, 76, 0, 34, 3, 17, 74, 90, 216, 185, 151, 76, 249, 232, 130, 154, 223, 42, 105, 40, 164, 79, 164, 47, 1, 140, 202, 90, 116, 124, 127, 97, 230, 207, 241, 144, 60, 131, 252, 38, 124, 83, 149, 31, 228, 131, 86, 134, 39, 201]))
dongle = getDongle(True)
publicKey = dongle.exchange(bytes("8004000000".decode('hex')))
print "publicKey " + str(publicKey).encode('hex')
try:
	offset = 0
	while offset != len(textToSign):
		if (len(textToSign) - offset) > 255:
			chunk = textToSign[offset : offset + 255] 
		else:
			chunk = textToSign[offset:]
		if (offset + len(chunk)) == len(textToSign):
			p1 = 0x80
		else:
			p1 = 0x00
		apdu = bytes("8002".decode('hex')) + chr(p1) + chr(0x00) + chr(len(chunk)) + bytes(chunk)
                start = time.time()
		signature = dongle.exchange(apdu)
                end = time.time()
		offset += len(chunk)  	
	print "signature " + str(signature).encode('hex')
        print "time to sign: " + str(end-start)
	publicKey = PublicKey(bytes(publicKey), raw=True)
	signature = publicKey.ecdsa_deserialize(bytes(signature))
	print "verified " + str(publicKey.ecdsa_verify(bytes(textToSign), signature))
except CommException as comm:
	if comm.sw == 0x6985:
		print "Aborted by user"
	else:
		print "Invalid status " + str(comm.sw) 

