import multiprocessing
import Bio
import Bio.pairwise2


h=">T0958"
Ts="MNKKSKQQEKLYNFIIAKSFQQPVGSTFTYGELRKKYNVVCSTNDQREVGRRFAYWIKYTPGLPFKIVGTKNGSLLYQKIGINPCNNSTPSKGGDC"
Tss="CCCCCHHHHHHHHHHHHHHHHCCCEEEECHHHHHHHCCCCCCHHHHHHHHHHHHHHHHHCCCCCEEEEEEECEEEEEEECCCCCCCCCCCCCCCCC"

# pdbidL = ["1BK6","1D8D","1D8E","1DCE","1E96","1EJL","1FPP","1FT1","1FT2","1G3J","1H2T","1H2U","1H2V","1H6K","1HE8","1I7W","1I7X","1IB1","1IBR","1IQ1","1JCQ","1JCR","1JCS","1JDH","1JPW","1K8K","1KPS","1KZO","1KZP","1LD7","1LD8","1LDJ","1LDK","1LSH","1LTX","1LUJ","1M1E","1M5N","1MZC","1N4P","1N4Q","1N4R","1N4S","1N52","1N54","1N94","1N95","1N9A","1NI1","1NL4","1O1R","1O1S","1O1T","1O5M","1OCC","1OCO","1OCR","1OCZ","1PJM","1PJN","1QBK","1QBQ","1QGK","1QGR","1R8Q","1R8S","1RE0","1S63","1S64","1S9D","1SA4","1SA5","1T08","1TH1","1TN6","1TN7","1TN8","1TNB","1TNO","1TNU","1TNY","1TNZ","1TU3","1TYQ","1U2V","1U6G","1UKL","1UW4","1V18","1V54","1V55","1WA5","1X79","1X81","1XQS","1Z2C","1Z5S","2BED","2BKU","2C0J","2C0L","2C1M","2C1T","2CFH","2DYR","2DYS","2EIJ","2EIK","2EIL","2EIM","2EIN","2F0Y","2GL7","2GRN","2GRO","2GRP","2GRQ","2GRR","2H4M","2H6F","2H6G","2H6H","2H6I","2HYE","2IE3","2IE4","2IEJ","2IO2","2IO3","2IY0","2J3R","2J3T","2J3W","2JDQ","2JKR","2KZT","2L0T","2NPP","2NYL","2NYM","2O98","2OCC","2P8Q","2P9I","2P9K","2P9L","2P9N","2P9P","2P9S","2P9U","2PF4","2PKG","2PQN","2PQR","2Q5D","2R2L","2VGL","2WTK","2Y69","2ZIR","2ZIS","2ZXW","3ABK","3ABL","3ABM","3AG1","3AG2","3AG3","3AG4","3ASN","3ASO","3AXJ","3AXY","3DPY","3DW8","3DXK","3DXM","3E30","3E32","3E33","3E34","3E37","3EA5","3EG5","3EU5","3EUV","3FEY","3FGA","3GB8","3GC3","3GNI","3IFQ","3J9T","3JB9","3K7V","3K7W","3KSL","3KSQ","3KXC","3LWW","3M50","3OBV","3OUW","3OUX","3PJA","3PZ4","3Q73","3Q7A","3Q7F","3QB5","3R9A","3RSE","3SBT","3SFX","3SFY","3TJ3","3TX7","3U88","3UIN","3UIO","3UIP","3UKR","3UKU","3UKZ","3UL0","3UL1","3ULE","3W5K","3WG7","3X2Q","3ZEF","3ZHP","4AP2","4APF","4BPL","4CR2","4FZA","4GTM","4GTO","4GTP","4GTQ","4GTR","4I43","4I5L","4I5N","4ILG","4ILH","4JD2","4JGH","4KBQ","4KP3","4KXK","4KYO","4L9P","4LNB","4LNG","4LWZ","4MBG","4N3Y","4N3Z","4P6Z","4PVZ","4R10","4R11","4UAD","4UAE","4UAF","4UQI","4YDE","4YDO","4YI0","5B1A","5B1B","5B3S","5D2M","5F5S","5F5T","5F5U","5F5V","5F74","5GAN","5GMK","5GXW","5IY5","5JCZ","5K9S","5N6N","5O9Z","5OWU","5W97","5WAU","5WQL","5X19","5X1B","5X1F","5XDQ"]
pdbidL = ["9ANT","1AHD","1SAN","2HOA","1HOM","1WI3","1P7I","2HDD","1P7J","1ENH","1DU0","3HDD","2HOS","2HOT","1HDD","2JWT","2P81","1ZTR","1JGG","1B8I","1FTZ","1IC8","1LFB","2LFB","1S7E","1MIJ","1XPX","1PUF","1B72","2CRA","1X2N","2E1O","2CUF","2ECC","1UHS","1ZQ3","1X58","1BW5","2CQX","1X2M","1K61","1AKH","1LE8","1YRN","1MNM","1APL","1F43","1MH3","1MH4","1IG7","1E3O","1GT0","1HF0","1OCT","1CQT","1O4X","1POG","1HDP","1OCP","2CUE","1FJL","1LFU","1DU6","1AU7","2LKX","1FTT","1X41","1NK3","1NK2","1VND","1QRY","1WH7","1WH5","3NAR","2DJO","2ECB","3A02","3LNQ","3A01","2R5Z","2R5Y","4UUT","3CMY","3RKQ","3K2A","4J19","2MSY","1WJH","2LMD","2LP0","2MW8","2DA6","2HI3","1GDT","1ZR4","2GM4","1ZR2","1RET","1RES","1JJ6","1IJW","1JKO","1HCR","1JKR","1JKP","1JJ8","1JKQ","2EZL","2EZK","2EZI","2EZH","1TC3","1U78","1UG2","1A5J","1GVD","1GV5","1GUU","1GV2","1H88","1MBJ","1MBK","1MBH","1MBE","1MBG","1MBF","1IDZ","1IDY","1MSE","1MSF","1H89","2CQQ","2CQR","1WGX","2CRG","2CU7","1XC5","2CJJ","1FEX","2IW5","2UXX","3ZN0","3ZMZ","2UXN","3ZMS","3ZMT","3ZN1","3ZMU","3ZMV","2XAG","2XAH","5L3D","2XAS","4UVB","2XAJ","4UV8","4UXN","2XAF","4UVA","4UV9","4UVC","2XAQ","1OFC","2CKX","2QHB","2AJE","1H8A","4A69","1W0T","1IV6","1ITY","1BA5","1W0U","1XG1","1VF9","1VFC","3SJM","1PDN","1K78","1MDM","6PAX","1IGN","3UKG","1IUF","1HLV","1BW6","1BL0","1D5Y","1UI5","1UI6","1T56","5NIO","3TP3","5NIM","5MYN","3TP0","5MYS","5MYT","5MXV","4M3D","5J3L","5MYL","5MYW","1U9N","5J1U","5J1Y","5MYR","4M3B","5J1R","5MXK","5NIZ","5MWO","4M3F","4M3G","5NJ0","4M3E","5MYM","1U9O","5IOY","5IPA","5IP6","5IOZ","2FX0","4JYK","1PB6","4XK4","3LOC","1RKT","1VI0","3BR0","1JT6","3BR3","3PM1","1RKW","2HQ5","3BTL","3BTI","2DTZ","3BT9","1JUS","1JTX","2G0E","3BTJ","1RPW","1QVT","1QVU","3BR5","3BTC","1JT0","2GBY","1JTY","1JUP","1JUM","2GEN","2FD5","2GFN","2D6Y","2OI8","2G3B","2G7S","2G7G","2HKU","2I10","2NP3","3C07","2HYJ","2ID3","2G7L","1SGM","1T33","2XPW","2XPV","2VKE","2XPU","2XPS","2XPT","4ABZ","4AUX","2X6O","1BJZ","2XB5","2O7O","2X9D","1ORK","1A6I","3FK7","1DU7","3FK6","2TCT","2FJ1","4V2F","1QPI","4V2G","1BJ0","2NS7","1BJY","2TRT","2IU5","2FQ4","1ZK8","2O7T","2YVE","1V7B","2ZOY","2DH0","2ZOZ","1Z0X","2FBQ","2NP5","2ID6","2IEK","1ZKG","3IH4","3IH2","3IH3","1Z77","4I76","2JJ7","2JK3","1IRZ","1NTC","4IHV","4IHW","4IHX","4IHY","3JR9","3JRH","3IV5","3JRG","3JRF","3JRC","3JRD","3JRI","3JRA","3JRE","3JRB","1ETX","1ETO","1ETV","1ETW","1ETY","1ETK","3FIS","1FIP","4FIS","1FIA","1F36","1ETQ","1UMQ","1G2H","1RR7","2COB","2GA1","2AO9","2DW4","2Z3Y","2Z5U","2EJR","2H94","2COM","2FQ3","2AQF","2CUJ","2AQE","2L3D","2JN6","2O3F","2NOG","2WV1","5GP9","5GPA","3WHB","3WHC","3ZOB","4B3A","4D7M","2XRL","5FKM","4B1R","2VKV","5FKK","4D7N","5FKO","5FKL","4AC0","5FKN","2ELK","2MG4","2H1K","3A03","4RDU","5GJK","2L7F","2L7M","2M0C","2LTP","2M34","2N8G","2K40","2DMU","2ME6","2DA2","1WQI","2DMQ","2DJN","2ME0","1YZ8","2DMT","2YUM","2DA3","2DA1","2DA5","2DN0","2L7Z","4D5F","2XSD","2VI6","2LD5","2DMS","2D9A","3Q0W","3Q0U","3Q0V","3O8G","3SDG","3O8H","3G1M","3Q3S","5EYR","4DW6","3G1L","3G1O","3SFI","5F1J","5F27","5EZH","5F0F","5F04","5EZG","5F0C","5F08","5F0H","2MGQ","2D5V","2VPR","4D5C","2OFL","3OSG","3OSF","2KDZ","2K9N","1SFE","1QNT","1EH7","1EH6","1EH8","1YFH","1T38","1T39","1MGT","3L00","3KZZ","3KZY","1C20","1KQQ","1IG6","2OEH","1RYU","1KKX","1KN5","2LI6","2CXY","2EH9","4LJX","2KK0","2LM1","2JRZ","2EQY","1BIA","1HXD","1BIB","2EWN","1J5Y","1JHF","1JHH","1LEA","1LEB","1B4A","2P5K","2P5L","1F9N","1AOY","5CJ9","1I5Z","3RYP","4FT8","1I6X","1HW5","4N9I","2CGP","3QOP","4N9H","1ZRF","1G6N","2GZW","3N4M","1RUN","1RUO","1O3R","1LB2","3RYR","1O3Q","1J59","1DBC","1O3S","1ZRE","1ZRD","1ZRC","1O3T","1CGP","2H6C","3E5U","3E6C","3E5X","3E6B","2H6B","3E6D","1FT9","2OZ6","1OMI","1ZYB","2GAU","2ZCW","2COH","2BGC","2BEO","1U2W","1R1U","1R1V","2M30","1Y0U","1R1T","1R23","1SMT","1R22","3F72","2KJC","2KJB","4GGG","1HW1","1E2X","1H9G","1HW2","1H9T","2HS5","3BWG","1V4R","1BM9","1F4K","1J0R","2EFW","2DPD","2DQR","2DPU","1B9M","1B9N","1O7L","1BJA","1I1S","1REP","2Z9O","2NRA","1HKQ","1FNN","1W5S","1W5T","1IN4","1J7K","1IN6","1IN7","1IN8","1IN5","1IXS","1HQC","1IXR","2FNA","2FOK","1FOK","1YQA","1UST","1USS","1UHM","1GHC","1HST","1D5V","3L2C","1E17","2A07","2HFH","2HDC","1KQ8","1JXS","3BPY","2C6Y","3CO6","1VTN","3COA","3QRF","3CO7","2K86","2D2W","2A3S","2BBY","1BBY","1DPU","1Z1D","4OU0","1CF7","1D8K","1D8J","1SFU","1J75","1OYI","1XMK","1QBJ","2GXB","3F21","3F22","3F23","2ACJ","1QGP","3IRR","3IRQ","2L54","2HEO","1DP7","1KA8","1WWX","1DUX","1GVJ","2STT","2STW","1MD0","1K79","1K7A","1R36","1ETC","1FLI","1AWC","1YO5","1BC8","1BC7","1K6O","1HBX","1PUE","2NNY","3WTS","3WTT","3WU1","3WU0","3WTZ","2YPR","5E8G","3WTU","3WTV","3WTY","3WTW","3RI4","1STW","3JTG","1HKS","1HKT","2HTS","3HTS","1FBQ","1FBS","1FYM","1FYL","1FBU","1FYK","3HSF","5HDN","5HDG","2LDU","1IF1","2PI0","1T2K","2O6G","2IRF","1IRF","1IRG","3QU6","2DTR","2QQA","2QQ9","1BI1","1BI0","1FWZ","1P92","1XCV","2QQB","1BI3","1BI2","2TDX","1G3T","1G3W","1G3S","1F5T","1G3Y","1DDN","1DPR","1C0W","2ISY","1FX7","2IT0","1U8R","2ISZ","1B1B","2EV0","1ON2","2EV6","1ON1","2F5D","2EV5","2F5E","2F5F","2F5C","2HYF","2HYG","4HX4","3R60","3R61","4HX7","4HV5","4HX8","4HV6","1B6A","1QZY","1BN5","1YW9","1B59","1BOA","1KQ9","1YW7","1R58","1KQ0","2GA2","1R5G","2OAZ","1R5H","2ADU","2EA4","2EA2","1YW8","1XGS","1WKM","1XGN","1XGM","1XGO","1T0F","1F1Z","1UB9","3MEX","3ECH","1LNW","1JGS","1Z91","1Z9C","2FRH","2FNP","1FZP","2HR3","2FBI","2FXA","2ETH","1S3J","1HSJ","1P4X","3CTA","2FBK","2BV6","3BRO","2FBH","1LJ9","3DEU","3QPT","3Q5F","2A61","5DD8","3VB2","3VOE","4JBA","3VOD","4ZZL","1FZN","4AIH","1QZZ","1R00","1XDS","1XDU","1KYZ","1KYW","5EEH","1TW3","5JR3","4WXH","1TW2","5EEG","1FP1","1FPQ","1FP2","1FPX","1I27","1J2X","1NHA","1ONV","2CSO","1W4M","1UHW","1V3F","1O7F","1FSH","5SUZ","5SUY","5LNP","1I1G","2ZNZ","2E1C","1RI7","2ZNY","2CG4","2CFX","1JMR","1MKM","1LDD","1LDJ","1U6G","4P5O","1LDK","1IUY","2HYE","1LVA","2UWM","1WSU","2PLY","2V9V","1KU9","1IXC","1IZ1","2ESN","2DT5","1XCB","1R72","2G9W","1SD7","1SAX","1OKR","1SD6","2D45","1P6R","2P7C","1SD4","1XSD","1P41","1P4A","1O57","1Q1H","1MZB","1OYW","1OYY","2P6R","2AXL","1PP7","1PP8","1UCR","1WQ2","2CQK","1ZH5","2VOO","2VON","2VOD","1YTY","2VOP","1S7A","1S29","4LCT","1UFM","3CHM","1RZ4","1WI9","4CR2","3TXN","3TXM","1TBX","2CO5","2PG4","1R7J","1XSX","1SFX","2D1H","1STZ","2ZXX","1WLQ","1XB4","1W7P","1U5T","1YLF","1XD7","1TQI","1TQM","1TQP","1ZAR","1ZAO","1W1W","1ULY","2CWE","1WJ5","1T6S","1YG2","2ESH","1XMA","1XN7","4AWX","2K02","3BP8","1Z6R","2HOE","1Z05","1XNP","2P4W","1T98","2B0L","2J5P","2J5O","2VE8","2VE9","5DEQ","5BS6","5DD4","2FML","1Z7U","2F2E","2FSW","2HZT","1YYV","2OD5","2OBP","2P8T","2HTJ","2IPQ","3BZ6","2NR3","2NS0","2HGC","2VQC","2CMX","3F0P","3F0O","3F2G","3F2F","5U7A","5U79","5U83","3F2H","5U7C","5U88","5C0U","5U82","3FN8","5C0T","5DSF","5U7B","1S6L","2A5Y","2DOA","2E5N","2GMG","1ZEL","2V7F","2DK5","2DK8","3NFI","3KX2","5JPT","2P6U","4ESF","3BJA","4ESB","2FE3","4A5N","2RGV","4A5M","3F8N","5HS7","5HS8","5HS9","5X14","5X13","5X12","5Y8T","5H20","3PFI","4ETS","2FMY","3REO","5CVV","5CVJ","5CVU","3TKY","3JW4","3GLX","3R6S","4BYY","4Q48","4Q47","4EJO","3RI2","5U8J","3HHH","1DB8","1DB7","1DB9","4I0B","3KCC","4HZF","4R8H","4I02","4I09","4I0A","3ROU","4I01","3JSO","3FWE","3RPQ","3RDI","3JSP","2WC2","4WF2","4MTE","4MTD","5VYV","2XIG","4CO8","5ILU","5HDK","5ILS","4IRG","4IRH","5ILV","4AVP","3CUQ","2AS5","2ZME","5OCN","2LNB","2DLL","2DAO","3K69","3F8B","3F8C","3F8F","4ZZD","5NL9","4E70","4EMS","4EVI","3P9C","3P9I","3P9K","4RAY","4RB0","4RAZ","4RB3","4RB1","4RB2","2QYO","1ZGJ","1ZGA","1ZG3","1ZHF","3BPV","3BPX","3QU3","2MD5","2LF7","3GW2","2KKO","3D0S","4A2U","3I54","3I59","2O03","2LKP","3H3U","3MZH","2JSC","2P5V","2P6S","2P6T","2MBF","5ERI","5DYM","6BLB","4GCV","2CYY","5FD5","5FD6","3CJN","3E6M","4OMY","4OMZ","2M5W","4PGG","4PGH","5FFX","4L9N","5F6F","4L9T","4LD5","4L9V","4HQM","5FFZ","4LLL","3ECO","4GXO","3HSR","3HRM","3HSE","4HBL","3KEO","3KEQ","3L7W","4LMY","4I7H","3EYY","3MWM","4HW0","2E7X","2E7W","2EFP","2PN6","2PMH","2YX4","2EFN","2YX7","2EFQ","2EFO","5ICE","5ICF","5ICC","5ICG","3ELK","3DF8","3K2Z","3IKV","4K2E","4OOI","2W57","3JTH","2FA5","3IWZ","3PQK","3PQJ","1OPC","1ODD","2JPB","1GXQ","1GXP","2Z33","1QQI","1KGS","2FF4","2FEZ","1P2F","1YS6","1YS7","1FSE","1A04","1JE8","1RNL","1ZG5","1ZG1","1L3L","1H0M","1YIO","1ZN2","1P4W","1FC3","1LQ1","2NAZ","3ULQ","2D1V","2KRF","3Q9V","3Q9S","2HWV","3ZQ7","4UHT","5JU7","2K4J","2HQR","2HQN","2M87","4NHJ","2MLK","2JZY","3C57","1ZLJ","1ZLK","2PMU","2OQR","2JPC","2ZXJ","2Z9M","4QWQ","2RNJ","4U88","4IXA","5DCM","3RJP","1HC8","1Y39","1QA6","1FOX","1FOW","1FOY","2FOW","1ACI","3CF5","2ZJQ","2ZJP","1XBP","2QAO","2QAM","2QBE","2QBG","2I2T","2I2V","2QOZ","2QP1","3DF4","3DF2","2QBC","2QBA","1VS8","1VS6","2AWB","2AW4","2QOX","2QOV","2QBI","2QBK","2VHM","2VHN","2Z4N","2Z4L","2RDO","3DEG","2J28","2GYC","2GYA","1VQ8","1VQO","1VQP","1YHQ","1S72","1VQM","1VQL","1VQK","1VQN","1YIJ","1YI2","1VQ9","1VQ5","1VQ4","1VQ6","1YIT","3CC2","1YJN","1YJ9","3CC7","1VQ7","3CCM","2QA4","3CC4","2OTL","3CCE","3CD6","3CCU","3CMA","3CCL","3CCV","3CCQ","2OTJ","2QEX","1YJW","3CCS","3CME","3CCR","3CCJ","1MMS","3CJR","3CJT","3CJQ","2NXN","2H8W","2E35","2E36","2E34","2HGQ","2HGJ","2HGU","1YL3","2B66","2B9N","2B9P","2QAN","2QAL","2QBD","2QBF","2QOY","2QP0","2I2P","3DF1","3DF3","2I2U","1VS5","1VS7","2AVY","2AW7","2QB9","2QBB","2VHO","2VHP","2QOU","2QOW","2QBH","2QBJ","2Z4K","2Z4M","1J5E","1FJG","2VQE","1XMQ","1N32","1XNQ","1XNR","1HR0","2UUB","1HNZ","1XMO","1I94","2UXD","2VQF","2UXC","2UUA","1HNW","1HNX","2E5L","2J02","2J00","2UUC","2UU9","1N33","2HHH","3D5C","3D5A","2UXB","1N34","1I96","2V46","2V48","2F4V","1N36","1I97","1I95","2QNH","2HGP","2HGI","2HGR","2OW8","1YL4","2B9O","1X18","2B64","2B9M","1G1X","1WHU","1E3P","1E3H","1K6Y","1WJF","1WJC","1WJB","1WJD","1WJE","1WJA","1E0E","1AUB","1TWF","1I50","1K83","1I3Q","1TWC","1TWA","4C2M","1TWH","1I6H","4C3H","1TWG","2NVY","2E2I","2NVT","2E2J","2VUM","2NVX","2B63","2YU9","2JA5","2JA7","2R92","2JA8","2R7Z","2R93","2JA6","2E2H","2NVZ","2B8K","3S14","1EF4","3CQZ","3M3Y","2NVQ","3S1N","3S1M","3RZO","4C3I","1JHG","3SSX","3SSW","2WRP","2OZ9","1P6Z","1TRO","3WRP","1WRP","1ZT9","1TRR","1MI7","1RCS","1WRS","1WRT","1CO0","2XDI","5TM0","1L8Q","2HCB","1J1V","2OA4","2JRT","1KU2","1RP3","1SC5","1SMY","2A6H","2BE5","2A68","2A69","1IW7","3EQL","2A6E","5TMC","2CW0","2P7V","4JK1","1TLH","1TTY","2K6X","1KU3","1RIO","1KU7","1OR7","2H27","1L0O","4G6D","4G94","4G8X","1XSV","1S7O","4G7H","4G7O","3DXJ","5TMF","1ZYR","4OIO","4Q4Z","4Q5S","1VZ0","1R71","1RQ6"]

subd = {('H','H'):3, ('E','E'):3, ('C','C'):1, ('H','E'):-3, ('E','H'): -3, ('H','C'):-2, ('E','C'):-2, ('C','H'): -2, ('C','E'): -2}

c = open("/home/tsukasa/casp_database/ss.txt").read()[1:].split("\n>")
def f(seq_c, ss_c):
    seq_head = seq_c.splitlines()[0]
    ss_head = ss_c.splitlines()[0]
    pdbid, chainid, _ = seq_head.split(":")
    if pdbid not in pdbidL:
        return None
    seq = "".join(seq_c.splitlines()[1:])
    ss = "".join(ss_c.splitlines()[1:]).translate(str.maketrans('HBEGITS ', 'HEEHHCCC'))
    score = Bio.pairwise2.align.localdd(Tss, ss, subd, -5,-2,-5,-2, score_only=True)
    return score, seq, ss, pdbid, chainid

with multiprocessing.Pool(processes=150) as pool:
    ret = pool.starmap(f, zip(c[::2],c[1::2]))
    ret = [e for e in ret if e is not None]

r = sorted(ret, reverse=True)[:100]

for rr in r:
    aligned_ss1,aligned_ss2,_,_,_ = Bio.pairwise2.align.localdd(Tss, rr[2], subd, -5,-2,-5,-2, one_alignment_only=True)[0]

    if len(aligned_ss1) > len(Tss):
        positions_where_addtional_gap_inserted_by_alignment = []
        index_on_ss = 0
        for i, aligned_ss1_r in enumerate(aligned_ss1):
            if len(Tss) <= index_on_ss:
                # end flanking gap tekina monowo soutei, koreizyou index_on_ss wo huyasenai joutai
                positions_where_addtional_gap_inserted_by_alignment.append(i)
                continue
            if aligned_ss1_r == Tss[index_on_ss]:
                index_on_ss += 1
                continue
            else:
                # start flanking gap tekina monowo soutei
                positions_where_addtional_gap_inserted_by_alignment.append(i)
        l_Ts = list(Ts)
        for i in  positions_where_addtional_gap_inserted_by_alignment:
            l_Ts.insert(i, "-")
        aTs = "".join(l_Ts)
    else:
        aTs = Ts

    if len(aligned_ss2) > len(rr[2]):
        positions_where_addtional_gap_inserted_by_alignment = []
        index_on_ss = 0
        for i, aligned_ss2_r in enumerate(aligned_ss2):
            if len(rr[2]) <= index_on_ss:
                # end flanking gap tekina monowo soutei, koreizyou index_on_ss wo huyasenai joutai
                positions_where_addtional_gap_inserted_by_alignment.append(i)
                continue
            if aligned_ss2_r == rr[2][index_on_ss]:
                index_on_ss += 1
                continue
            else:
                # start flanking gap tekina monowo soutei
                positions_where_addtional_gap_inserted_by_alignment.append(i)
        l_ts = list(rr[1])
        for i in  positions_where_addtional_gap_inserted_by_alignment:
            l_ts.insert(i, "-")
        a_ts = "".join(l_ts)
    else:
        a_ts  = rr[1]


    print(rr[3], rr[4], rr[0])
    # query_id = f"T0958_N4C10-ss_l_3-{rr[3]}_{rr[4]}"
    # print(aTs)
    # print(aligned_ss1)
    # print(aligned_ss2)
    # print(a_ts)

#     pir = f"""
# >P1;{query_id}
# sequence:{query_id}: : : : : : :0.00:0.00
# {aTs}
# *
# >P1;{rr[3]}_{rr[4]}
# structure:{rr[3].lower()}:.:{rr[4]}:.:{rr[4]}: : :-1.00:-1.00
# {a_ts}
# *
# """
#     pir_path="pir2/"+query_id+".pir"
#     with open(pir_path, "w") as f:
#         f.write(pir)