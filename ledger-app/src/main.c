/*******************************************************************************
*   Ledger Blue
*   (c) 2016 Ledger
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*      http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
********************************************************************************/

#include "os.h"
#include "cx.h"
#include "bls-embedded.h"

#include "os_io_seproxyhal.h"

unsigned char G_io_seproxyhal_spi_buffer[IO_SEPROXYHAL_BUFFER_SIZE_B];

static unsigned int current_text_pos; // parsing cursor in the text to display
static unsigned int text_y;           // current location of the displayed text

// UI currently displayed
enum UI_STATE { UI_IDLE, UI_TEXT, UI_APPROVAL };

enum UI_STATE uiState;

// Specific to Ledger nano X
#include "ux.h"
ux_state_t G_ux;
bolos_ux_params_t G_ux_params;

static const bagl_element_t *io_seproxyhal_touch_exit(const bagl_element_t *e);
static const bagl_element_t*
io_seproxyhal_touch_approve(const bagl_element_t *e);
static const bagl_element_t *io_seproxyhal_touch_deny(const bagl_element_t *e);

static void ui_idle(void);
//static unsigned char display_text_part(void);
static void ui_text(void);
static void ui_approval(void);

#define MAX_CHARS_PER_LINE 49
#define DEFAULT_FONT BAGL_FONT_OPEN_SANS_LIGHT_16px | BAGL_FONT_ALIGNMENT_LEFT
#define SIG_BYTES 96
#define PUBKEY_BYTES 192

#define CLA 0xe0
#define INS_SIGN 0x02
#define INS_GET_PUBLIC_KEY 0x04
#define INS_GET_APP_CONFIGURATION 0x06
#define P1_LAST 0x80
#define P1_MORE 0x00

#define APP_TYPE 0x04
#define LEDGER_MAJOR_VERSION 1
#define LEDGER_MINOR_VERSION 0
#define LEDGER_PATCH_VERSION 0

//static char lineBuffer[MAX_CHARS_PER_LINE+1];

UX_STEP_NOCB(
    ux_idle_flow_1_step,
    pnn,
    {
      &C_icon_certificate,
      "Waiting for message",
    });
UX_STEP_VALID(
    ux_idle_flow_2_step,
    pbb,
    io_seproxyhal_touch_exit(NULL),
    {
      &C_icon_crossmark,
      "Exit app",
    });

UX_FLOW(ux_idle_flow,
  &ux_idle_flow_1_step,
  &ux_idle_flow_2_step
);

///////////////////////////////////

UX_STEP_NOCB(
    ux_approve_flow_1_step,
    pnn,
    {
      &C_icon_certificate,
      "Waiting for message",
    });
UX_STEP_VALID(
    ux_approve_flow_2_step,
    pbb,
    io_seproxyhal_touch_approve(NULL),
    {
      &C_icon_validate_14,
      "Sign message",
    });
UX_STEP_VALID(
    ux_approve_flow_3_step,
    pbb,
    io_seproxyhal_touch_deny(NULL),
    {
      &C_icon_crossmark,
      "Deny signature",
    });

UX_FLOW(ux_approval_flow,
  &ux_approve_flow_1_step,
  &ux_approve_flow_2_step,
  &ux_approve_flow_3_step
);


// Text flow displays message for debugging
////////////////////////////////////
UX_STEP_NOCB(
    ux_text_flow_1_step,
    pnn, 
    {
      &C_icon_certificate,
      "This is the text flow",
    });

UX_STEP_NOCB(
    ux_text_flow_2_step,
    bnnn_paging,
    {
      .title = "Message",
//      .text = lineBuffer,
    });    

UX_STEP_VALID(
    ux_text_flow_3_step,
    pbb,
    io_seproxyhal_touch_approve(NULL),
    {
      &C_icon_validate_14,
      "Sign message",
    });

UX_STEP_VALID(
    ux_text_flow_4_step,
    pbb,
    io_seproxyhal_touch_deny(NULL),
    {
      &C_icon_crossmark,
      "Reject message",
    });

UX_FLOW(ux_text_flow,
  &ux_text_flow_1_step,
  &ux_text_flow_2_step,
  &ux_text_flow_3_step,
  &ux_text_flow_4_step
);

static const bagl_element_t *io_seproxyhal_touch_exit(const bagl_element_t *e) {
    // Go back to the dashboard
    os_sched_exit(0);
    return NULL; // do not redraw the widget
}

static void gen_private_key(uint64_t *out_private_key) {
    unsigned int bip32_path[5];
    bip32_path[0] = 44     | 0x80000000; // BIP44 specification
    bip32_path[1] = 0x8000ce10; // Specifies Celo. Note: BIP44 path not publicly secured yet
    bip32_path[2] = 0      | 0x80000000; // Account value
    bip32_path[3] = 0;
    bip32_path[4] = 0;                   // Index of derived child key 
    os_perso_derive_node_bip32(CX_CURVE_256K1, bip32_path,
                               sizeof(bip32_path) / sizeof(bip32_path[0]), 
                               (unsigned char*)out_private_key,
                               NULL);
    // The modulus of the scalar field of 256k1 is larger than for bls12_377. If the entropy generated does not fit
    // into an element of the bls12_377 scalar field, we increment the index in the bip32 path and try again
    if (!is_valid_key((const uint8_t*)out_private_key)) {
        bip32_path[4] += 1;
        os_perso_derive_node_bip32(CX_CURVE_256K1, bip32_path,
                                   sizeof(bip32_path) / sizeof(bip32_path[0]),
                                   (unsigned char*)out_private_key,
                                   NULL);
    }
}

static const bagl_element_t*
io_seproxyhal_touch_approve(const bagl_element_t *e) {
    uint64_t in_private_key[4];
    gen_private_key(in_private_key);
    uint8_t *input_hash = G_io_apdu_buffer + 5;
    sign_hash(in_private_key, input_hash, G_io_apdu_buffer);
    G_io_apdu_buffer[SIG_BYTES] = 0x90;
    G_io_apdu_buffer[SIG_BYTES+1] = 0x00;
    // Send back the response, do not restart the event loop
    io_exchange(CHANNEL_APDU | IO_RETURN_AFTER_TX, SIG_BYTES+2);
    // Display back the original UX
    ui_idle();
    return 0; // do not redraw the widget
}

static const bagl_element_t *io_seproxyhal_touch_deny(const bagl_element_t *e) {
    G_io_apdu_buffer[0] = 0x69;
    G_io_apdu_buffer[1] = 0x85;
    // Send back the response, do not restart the event loop
    io_exchange(CHANNEL_APDU | IO_RETURN_AFTER_TX, 2);
    // Display back the original UX
    ui_idle();
    return 0; // do not redraw the widget
}

unsigned short io_exchange_al(unsigned char channel, unsigned short tx_len) {
    switch (channel & ~(IO_FLAGS)) {
    case CHANNEL_KEYBOARD:
        break;

    // multiplexed io exchange over a SPI channel and TLV encapsulated protocol
    case CHANNEL_SPI:
        if (tx_len) {
            io_seproxyhal_spi_send(G_io_apdu_buffer, tx_len);

            if (channel & IO_RESET_AFTER_REPLIED) {
                reset();
            }
            return 0; // nothing received from the master so far (it's a tx
                      // transaction)
        } else {
            return io_seproxyhal_spi_recv(G_io_apdu_buffer,
                                          sizeof(G_io_apdu_buffer), 0);
        }

    default:
        THROW(INVALID_PARAMETER);
    }
    return 0;
}

static void sample_main(void) {
    volatile unsigned int rx = 0;
    volatile unsigned int tx = 0;
    volatile unsigned int flags = 0;


    // next timer callback in 500 ms
    UX_CALLBACK_SET_INTERVAL(500);

    // DESIGN NOTE: the bootloader ignores the way APDU are fetched. The only
    // goal is to retrieve APDU.
    // When APDU are to be fetched from multiple IOs, like NFC+USB+BLE, make
    // sure the io_event is called with a
    // switch event, before the apdu is replied to the bootloader. This avoid
    // APDU injection faults.
    for (;;) {
        volatile unsigned short sw = 0;

        BEGIN_TRY {
            TRY {
                rx = tx;
                tx = 0; // ensure no race in catch_other if io_exchange throws
                        // an error
                rx = io_exchange(CHANNEL_APDU | flags, rx);
                flags = 0;

                // no apdu received, well, reset the session, and reset the
                // bootloader configuration
                if (rx == 0) {
                    THROW(0x6982);
                }

                if (G_io_apdu_buffer[0] != CLA) {
                    THROW(0x6E00);
                }

                switch (G_io_apdu_buffer[1]) {
                case INS_SIGN: {
                    if (G_io_apdu_buffer[2] != P1_LAST) {
                        THROW(0x6A86);
                    }
                    // Wait for the UI to be completed
                    current_text_pos = 0;
                    text_y = 60;
                    G_io_apdu_buffer[5 + G_io_apdu_buffer[4]] = '\0';

                    //display_text_part();
                    ui_text();

                    flags |= IO_ASYNCH_REPLY;
                } break;

                case INS_GET_PUBLIC_KEY: {
	            uint64_t private_key[4];
                    gen_private_key(private_key);					 
		    uint8_t public_key[PUBKEY_BYTES];
		    get_pubkey(private_key, public_key);
		    os_memmove(G_io_apdu_buffer, public_key, PUBKEY_BYTES);
		    tx = PUBKEY_BYTES;
		    THROW(0x9000);
                } break;

		case INS_GET_APP_CONFIGURATION: {
		    G_io_apdu_buffer[0] = 0x00;
		    G_io_apdu_buffer[1] = APP_TYPE;
		    G_io_apdu_buffer[2] = LEDGER_MAJOR_VERSION;
		    G_io_apdu_buffer[3] = LEDGER_MINOR_VERSION;
		    G_io_apdu_buffer[4] = LEDGER_PATCH_VERSION;
		    tx = 5;
		    THROW(0x9000);
		} break;

                case 0xFF: // return to dashboard
                    goto return_to_dashboard;

                default:
                    THROW(0x6D00);
                    break;
                }
            }
            CATCH_OTHER(e) {
                switch (e & 0xF000) {
                case 0x6000:
                case 0x9000:
                    sw = e;
                    break;
                default:
                    sw = 0x6800 | (e & 0x7FF);
                    break;
                }
                // Unexpected exception => report
                G_io_apdu_buffer[tx] = sw >> 8;
                G_io_apdu_buffer[tx + 1] = sw;
                tx += 2;
            }
            FINALLY {
            }
        }
        END_TRY;
    }

return_to_dashboard:
    return;
}

void io_seproxyhal_display(const bagl_element_t *element) {
    io_seproxyhal_display_default((bagl_element_t *)element);
}

// Pick the text elements to display
/*static unsigned char display_text_part() {
    unsigned int i;
    WIDE char *text = (char*) G_io_apdu_buffer + 5;
    if (text[current_text_pos] == '\0') {
        return 0;
    }
    i = 0;
    while ((text[current_text_pos] != 0) && (text[current_text_pos] != '\n') &&
           (i < MAX_CHARS_PER_LINE)) {
        lineBuffer[i++] = text[current_text_pos];
        current_text_pos++;
    }
    if (text[current_text_pos] == '\n') {
        current_text_pos++;
    }
    lineBuffer[i] = '\0';
    return 1;
}*/

static void ui_idle(void) {
    uiState = UI_IDLE;
    ux_flow_init(0, ux_idle_flow, NULL);
}

static void ui_text(void) {
    uiState = UI_TEXT;
    // Display message if debug enabled
    #ifdef DEBUG
        ux_flow_init(0, ux_text_flow, NULL);
    #else
        io_seproxyhal_touch_approve(NULL);
    #endif
}

static void ui_approval(void) {
    uiState = UI_APPROVAL;
    ux_flow_init(0, ux_approval_flow, NULL);
}

unsigned char io_event(unsigned char channel) {
    // nothing done with the event, throw an error on the transport layer if
    // needed

    // can't have more than one tag in the reply, not supported yet.
    switch (G_io_seproxyhal_spi_buffer[0]) {
    case SEPROXYHAL_TAG_FINGER_EVENT:
        UX_FINGER_EVENT(G_io_seproxyhal_spi_buffer);
        break;

    case SEPROXYHAL_TAG_BUTTON_PUSH_EVENT:
        UX_BUTTON_PUSH_EVENT(G_io_seproxyhal_spi_buffer);
        break;

    case SEPROXYHAL_TAG_DISPLAY_PROCESSED_EVENT:
        if ((uiState == UI_TEXT) &&
            (os_seph_features() &
             SEPROXYHAL_TAG_SESSION_START_EVENT_FEATURE_SCREEN_BIG)) {
        //    if (!display_text_part()) {
        //        ui_approval();
        //    } else {
                UX_REDISPLAY();
        //    }
        } else {
            UX_DISPLAYED_EVENT();
        }
        break;

    // unknown events are acknowledged
    default:
        UX_DEFAULT_EVENT();
        break;
    }

    // close the event if not done previously (by a display or whatever)
    if (!io_seproxyhal_spi_is_status_sent()) {
        io_seproxyhal_general_status();
    }

    // command has been processed, DO NOT reset the current APDU transport
    return 1;
}

__attribute__((section(".boot"))) int main(void) {
    // exit critical section
    __asm volatile("cpsie i");

    current_text_pos = 0;
    text_y = 60;
    //hashTainted = 1;
    uiState = UI_IDLE;

    // ensure exception will work as planned
    os_boot();

    UX_INIT();
    if (G_ux.stack_count == 0) {
        ux_stack_push();
    }	

    BEGIN_TRY {
        TRY {
            io_seproxyhal_init();

#ifdef LISTEN_BLE
            if (os_seph_features() &
                SEPROXYHAL_TAG_SESSION_START_EVENT_FEATURE_BLE) {
                BLE_power(0, NULL);
                // restart IOs
                BLE_power(1, NULL);
            }
#endif

            USB_power(0);
            USB_power(1);

            ui_idle();

            sample_main();
        }
        CATCH_OTHER(e) {
        }
        FINALLY {
        }
    }
    END_TRY;
}
