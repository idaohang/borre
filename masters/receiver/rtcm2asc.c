/***************************************************************************
**                                                                        **
**                            RTCM2ASC.C                                  **
**                                                                        **
**----------------------------------------------------------------------- **
**                                                                        **
** Copyright (c) 2000 by Mikkel Jensen                                    **
**              Danish GPS Center                                         **
**              Aalborg University                                        **
**              Niels Jernes Vej 12                                       **
**              DK-9920 Aalborg East                                      **
**              Denmark                                                   **
**              email: mikkel@kom.auc.dk                                  **
**                                                                        **
** Permission is granted to copy and distribute this file in modified     **
** or unmodified form, for noncommercial use, provided (a) this copyright **
** notice is preserved, (b) no attempt is made to restrict redistribution **
** of this file, and (c) this file is not distributed as part of any      **
** collection whose redistribution is restricted by a compilation         **
** copyright.                                                             **
**                                                                        **
** This software is distributed as it is in the hope that it will be      **
** useful, but WITHOUT ANY WARRANTY; without even the implied warranty of **
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE                    **
**                                                                        **
**------------------------------------------------------------------------**
**                                                                        **
** history:                                                               **
**   rev    date    author  comment                                       **
**                                                                        **
**   1.0   001220   MJ      First Release, compiled with Visual C++ 6     **
**                                                                        **
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <conio.h>

int bitparity[6];
int bitword[24];
char rawbuf[255];
char ch;
int k;
signed char d;
int l;
int d29star = 0;
int d30star = 1;

int bool_cslip = 0;
int bool_dprc = 0;
int bool_drrc = 0;
int bool_factor = 0;
int bool_freq = 0;
int bool_gnsstime = 0;
int bool_gpsglo = 0;
int bool_health = 0;
int bool_iod = 0;
int bool_length = 0;
int bool_mpe = 0;
int bool_msg = 0;
int bool_mtype = 0;
int bool_multi = 0;
int bool_pcind = 0;
int bool_phase = 0;
int bool_pr = 0;
int bool_prc = 0;
int bool_pre = 0;
int bool_prn = 0;
int bool_qual = 0;
int bool_refid = 0;
int bool_rrc = 0;
int bool_seqno = 0;
int bool_smooth = 0;
int bool_udre = 0;
int bool_xyz = 0;
int bool_zcount = 0;

int count;

div_t div_result;

FILE *rtcmbin;
FILE *mtype1, *mtype2, *mtype3, *mtype6, *mtype16, *mtype18, *mtype19, *seqno;

/*****************************************************************************
 * BITPOINT                                                                  *
 *                                                                           *
 * This small function accepts a byte "bt" and an integer "bp" as its input  *
 * variables and returns the value of bit number "bp" in "bt" as an integer  *
 * either 0 or 1                                                             *
 *****************************************************************************/

int bitpoint(char bt, int bp){
	return (bt >> bp) & 0x01;
}

/*****************************************************************************
 * LOADWORD                                                                  *
 *                                                                           *
 * LOADWORD takes five bytes ch1-ch5 as its input variables. From these it   *
 * generates the corresponding GPS word taking into account both parity and  *
 * the Least Significant Bit First Rule. It uses BITPOINT to perform the     *
 * byte roll part. The result is written to the Bitword variable and the     *
 * has no return value                                                       *
 *****************************************************************************/

void loadword(char ch1, char ch2, char ch3, char ch4, char ch5){
 
	int cnt;

	if (d30star == 0){
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt] = bitpoint(ch1, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 6] = bitpoint(ch2, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 12] = bitpoint(ch3, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 18] = bitpoint(ch4, cnt);
		}
	}
	else {
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt] = 1 - bitpoint(ch1, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 6] = 1 - bitpoint(ch2, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 12] = 1 - bitpoint(ch3, cnt);
		}
		for (cnt = 0; cnt < 6; cnt++){
			bitword[cnt + 18] = 1 - bitpoint(ch4, cnt);
		}
	}
	for (cnt = 0; cnt < 6; cnt++){
		bitparity[cnt] = bitpoint(ch5, cnt);
	}

	d29star = bitparity[4];
	d30star = bitparity[5];
}

/*****************************************************************************
 * READMESSAGE1                                                              *
 *                                                                           *
 * This function reads and decodes type 1 message and writes the result to a *
 * text file. The input to the function is the length of the message in      *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/

void readmessage1(int msglen){
	int type1packet[40];
	int m;
	int n = 10;
	int o;
	int q = 0;
	int p = 0;
	unsigned char scale_factor = 0;
	unsigned char udre = 0;
	unsigned char sat_id = 0;
	signed short prc = 0;
	signed char rrc = 0;
	unsigned char iod = 0;
	int j;

	double sfprc;
	double sfrrc;

	for (m = 0; m < msglen; m++){

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		p = 0;

		for (o = 0; o < 3; o++){
			
			for (j = 0; j < 8; j++){
				type1packet[q] = bitword[p];
				q++;
				p++;
			}

			if (q == 40){  // packet done, write to file
				
				sat_id = 0;
				prc = 0;
				rrc = 0;
				iod = 0;

				scale_factor = type1packet[0];
				udre = type1packet[1]*2 + type1packet[2];
				for (j = 0; j < 5; j++){
					sat_id = sat_id | (type1packet[7 - j] << j);
				}
				for (j = 0; j < 16; j++){
					prc = prc | (type1packet[23 - j] << j);
				}
				for (j = 0; j < 8; j++){
					rrc = rrc | (type1packet[31 - j] << j);
				}
				for (j = 0; j < 8; j++){
					iod = iod | (type1packet[39 - j] << j);
				}

				if (scale_factor == 0){
					sfprc = 0.02;
					sfrrc = 0.002;
				}
				else {
					sfprc = 0.32;
					sfrrc = 0.032;
				}

				if (sat_id == 0){
					sat_id = 32;
				}

				if (bool_factor){
					fprintf(mtype1,"%d ",scale_factor);
				}

				if (bool_udre){
					fprintf(mtype1,"%d ",udre);
				}

				if (bool_prn){
					fprintf(mtype1,"%2.0f ",(double)sat_id);
				}

				if (bool_prc){
					fprintf(mtype1,"%9.2f ",prc*sfprc);
				}

				if (bool_rrc){
					fprintf(mtype1,"%6.3f ",rrc*sfrrc);
				}

				if (bool_iod){
					fprintf(mtype1,"%3.0f ",(double)iod);
				}
				
				//fprintf(mtype1,"\t%d\t%d\t%d\t%d\t%d\t%d",scale_factor,udre,sat_id,prc,rrc,iod);

				q = 0;
			}


		}
		
	
	}

	fprintf(mtype1,"\n");
	
}

/*****************************************************************************
 * READMESSAGE2                                                              *
 *                                                                           *
 * This function reads and decodes type 2 message and writes the result to a *
 * text file. The input to the function is the length of the message in      *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/

void readmessage2(int msglen){
	int type2packet[40];
	int m;
	int n = 10;
	int o;
	int q = 0;
	int p = 0;
	unsigned char scale_factor = 0;
	unsigned char udre = 0;
	unsigned char sat_id = 0;
	signed short dprc = 0;
	signed char drrc = 0;
	unsigned char iod = 0;
	int j;

	double sfdprc;
	double sfdrrc;

	for (m = 0; m < msglen; m++){

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		p = 0;

		for (o = 0; o < 3; o++){
			
			for (j = 0; j < 8; j++){
				type2packet[q] = bitword[p];
				q++;
				p++;
			}

			if (q == 40){  // packet done, write to file
				
				sat_id = 0;
				dprc = 0;
				drrc = 0;
				iod = 0;

				scale_factor = type2packet[0];
				udre = type2packet[1]*2 + type2packet[2];
				for (j = 0; j < 5; j++){
					sat_id = sat_id | (type2packet[7 - j] << j);
				}
				for (j = 0; j < 16; j++){
					dprc = dprc | (type2packet[23 - j] << j);
				}
				for (j = 0; j < 8; j++){
					drrc = drrc | (type2packet[31 - j] << j);
				}
				for (j = 0; j < 8; j++){
					iod = iod | (type2packet[39 - j] << j);
				}

				if (scale_factor == 0){
					sfdprc = 0.02;
					sfdrrc = 0.002;
				}
				else {
					sfdprc = 0.32;
					sfdrrc = 0.032;
				}

				if (sat_id == 0){
					sat_id = 32;
				}

				if (bool_factor){
					fprintf(mtype2,"%d ",scale_factor);
				}

				if (bool_udre){
					fprintf(mtype2,"%d ",udre);
				}

				if (bool_prn){
					fprintf(mtype2,"%2.0f ",(double)sat_id);
				}

				if (bool_dprc){
					fprintf(mtype2,"%9.2f ",dprc*sfdprc);
				}

				if (bool_drrc){
					fprintf(mtype2,"%6.3f ",drrc*sfdrrc);
				}

				if (bool_iod){
					fprintf(mtype2,"%3.0f ",(double)iod);
				}

				q = 0;
			}


		}
		
	
	}

	fprintf(mtype2,"\n");
	
}

/*****************************************************************************
 * READMESSAGE3                                                              *
 *                                                                           *
 * This function reads and decodes type 3 message and writes the result to a *
 * text file. The input to the function is the length of the message in      *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/


void readmessage3(int msglen){
	int type3packet[96];
	int m;
	int n = 10;
	int q = 0;
	int p;
	signed int xcoord;
	signed int ycoord;
	signed int zcoord;
	int j;

	for (m = 0; m < msglen; m++){
		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		for (p = 0; p < 24; p++){
			type3packet[q] = bitword[p];
			q++;
		}


	}

	xcoord = 0;
	ycoord = 0;
	zcoord = 0;

	for (j = 0; j < 32; j++){
		xcoord = xcoord | (type3packet[31 - j] << j);
	}
	
	for (j = 0; j < 32; j++){
		ycoord = ycoord | (type3packet[63 - j] << j);
	}

	for (j = 0; j < 32; j++){
		zcoord = zcoord | (type3packet[95 - j] << j);
	}

	if (bool_xyz){
		fprintf(mtype3,"%11.2f %11.2f %11.2f ",xcoord*0.01,ycoord*0.01,zcoord*0.01);
	}

	fprintf(mtype3,"\n");

}

/*****************************************************************************
 * READMESSAGE16                                                             *
 *                                                                           *
 * This function reads and decodes type 16 message and writes the result to  *
 * a text file. The input to the function is the length of the message in    *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/

void readmessage16(int msglen){
	char strmessage[90];
	char chr;
	int m;
	int n = 10;
	int q = 0;
	int p;
	int j;
	

	strmessage[0] = '\0';
	for (m = 0; m < msglen; m++){
		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;
		for (p = 1; p < 4; p++){
			chr = 0;
			for (j = 0; j < 8; j++){
				chr = chr | (bitword[(p*8 - 1) - j] << j);
			}
			strmessage[q] = chr;
			q++;
		}

	}

	strmessage[q] = '\0';

	if (bool_msg){
		fprintf(mtype16,"%s",strmessage);
	}

	fprintf(mtype16,"\n");

}

/*****************************************************************************
 * READMESSAGE18                                                             *
 *                                                                           *
 * This function reads and decodes type 18 message and writes the result to  *
 * a text file. The input to the function is the length of the message in    *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/

void readmessage18(int msglen){
	const double phasefactor = 0.00390625;
	int extheader[24];
	int type18packet[48];
	int n = 10;
	int j;
	int q = 0;
	int p;
	int m;
	unsigned char freq_indicator = 0;
	unsigned int gnss_time = 0;

	unsigned char multi_message;
	unsigned char pc_indicator;
	unsigned char gps_glonass;
	unsigned char sat_id;
	unsigned char data_qual;
	unsigned char cum_loss_indicator;
	signed int carrier_phase;

	loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
	n = n + 5;

	for (p = 0; p < 24; p++){
		extheader[p] = bitword[p];
	}

	for (j = 0; j < 2; j++){
		freq_indicator = freq_indicator | (extheader[1 - j] << j);
	}

	for (j = 0; j < 20; j++){
		gnss_time = gnss_time | (extheader[23 - j] << j);
	}

	if (bool_freq){
		fprintf(mtype18,"%d ",freq_indicator);
	}

	if (bool_gnsstime){
		fprintf(mtype18,"%6.0f ",(double)gnss_time);
	}

	for (m = 0; m < (msglen - 1)/2; m++){
		q = 0;
		sat_id = 0;
		data_qual = 0;
		cum_loss_indicator = 0;
		carrier_phase = 0;

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		for (p = 0; p < 24; p++){
			type18packet[q] = bitword[p];
			q++;
		}		

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		for (p = 0; p < 24; p++){
			type18packet[q] = bitword[p];
			q++;
		}

		multi_message = type18packet[0];
		pc_indicator = type18packet[1];
		gps_glonass = type18packet[2];

		for (j = 0; j < 5; j++){
			sat_id = sat_id | (type18packet[7 - j] << j);
		}

		if (sat_id == 0){
			sat_id = 32;
		}
		
		for (j = 0; j < 3; j++){
			data_qual = data_qual | (type18packet[10 - j] << j);
		}

		for (j = 0; j < 5; j++){
			cum_loss_indicator = cum_loss_indicator | (type18packet[15 - j] << j);
		}

		for (j = 0; j < 32; j++){
			carrier_phase = carrier_phase | (type18packet[47 - j] << j);
		}

		if (bool_multi){
			fprintf(mtype18,"%d ",multi_message);
		}
		
		if (bool_pcind){
			fprintf(mtype18,"%d ",pc_indicator);
		}

		if (bool_gpsglo){
			fprintf(mtype18,"%d ",gps_glonass);
		}

		if (bool_prn){
			fprintf(mtype18,"%2.0f ",(double)sat_id);
		}

		if (bool_qual){
			fprintf(mtype18,"%d ",data_qual);
		}

		if (bool_cslip){
			fprintf(mtype18,"%2.0f ",(double)cum_loss_indicator);
		}

		if (bool_phase){
			fprintf(mtype18,"%12.3f ",phasefactor*carrier_phase);
		}
		
	}

	fprintf(mtype18,"\n");

}

/*****************************************************************************
 * READMESSAGE19                                                             *
 *                                                                           *
 * This function reads and decodes type 19 message and writes the result to  *
 * a text file. The input to the function is the length of the message in    *
 * bytes. LOADWORD is used to read the GPS words. There is no return value   *
 *****************************************************************************/

void readmessage19(int msglen){
	const double prfactor = 0.02;
	int extheader[24];
	int type19packet[48];
	int n = 10;
	int j;
	int q = 0;
	int p;
	int m;
	unsigned char freq_indicator = 0;
	unsigned char smoothing = 0;
	unsigned int gnss_time = 0;

	unsigned char multi_message;
	unsigned char pc_indicator;
	unsigned char gps_glonass;
	unsigned char sat_id;
	unsigned char data_qual;
	unsigned char mp_error;
	unsigned int pseudo_range;

	loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
	n = n + 5;

	for (p = 0; p < 24; p++){
		extheader[p] = bitword[p];
	}

	for (j = 0; j < 2; j++){
		freq_indicator = freq_indicator | (extheader[1 - j] << j);
	}

	for (j = 0; j < 2; j++){
		smoothing = smoothing | (extheader[3 - j] << j);
	}

	for (j = 0; j < 20; j++){
		gnss_time = gnss_time | (extheader[23 - j] << j);
	}

	if (bool_freq){
		fprintf(mtype19,"%d ",freq_indicator);
	}

	if (bool_smooth){
		fprintf(mtype19,"%d ",smoothing);
	}

	if (bool_gnsstime){
		fprintf(mtype19,"%6.0f ",(double)gnss_time);
	}

	for (m = 0; m < (msglen - 1)/2; m++){
		q = 0;
		sat_id = 0;
		data_qual = 0;
		mp_error = 0;
		pseudo_range = 0;

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		for (p = 0; p < 24; p++){
			type19packet[q] = bitword[p];
			q++;
		}		

		loadword(rawbuf[n],rawbuf[n+1],rawbuf[n+2],rawbuf[n+3],rawbuf[n+4]);
		n = n + 5;

		for (p = 0; p < 24; p++){
			type19packet[q] = bitword[p];
			q++;
		}

		multi_message = type19packet[0];
		pc_indicator = type19packet[1];
		gps_glonass = type19packet[2];

		for (j = 0; j < 5; j++){
			sat_id = sat_id | (type19packet[7 - j] << j);
		}
		
		if (sat_id == 0){
			sat_id = 32;
		}

		for (j = 0; j < 4; j++){
			data_qual = data_qual | (type19packet[11 - j] << j);
		}

		for (j = 0; j < 4; j++){
			mp_error = mp_error | (type19packet[15 - j] << j);
		}

		for (j = 0; j < 32; j++){
			pseudo_range = pseudo_range | (type19packet[47 - j] << j);
		}

		if (bool_multi){
			fprintf(mtype19,"%d ",multi_message);
		}
		
		if (bool_pcind){
			fprintf(mtype19,"%d ",pc_indicator);
		}

		if (bool_gpsglo){
			fprintf(mtype19,"%d ",gps_glonass);
		}

		if (bool_prn){
			fprintf(mtype19,"%2.0f ",(double)sat_id);
		}

		if (bool_qual){
			fprintf(mtype19,"%2.0f ",(double)data_qual);
		}

		if (bool_mpe){
			fprintf(mtype19,"%2.0f ",(double)mp_error);
		}

		if (bool_pr){
			fprintf(mtype19,"%11.2f ",prfactor*pseudo_range);
		}
		
	}

	fprintf(mtype19,"\n");

}

/*****************************************************************************
 * READHEADER                                                                *
 *                                                                           *
 * This function reads and decodes the two word header information. The      *
 * output is written to file and from the informtion in the header the       *
 * correct function is called to decode the rest of the message. This        *
 * function also uses LOADWORD to read the GPS word                          *
 *****************************************************************************/


void readheader(void){
	int j;

	int preamble = 0;
	int message_type = 0;
	int station_id = 0;
	int mod_zcount = 0;
	int seq_no = 0;
	int length = 0;
	int health = 0;

	FILE *header;

	loadword(rawbuf[0],rawbuf[1],rawbuf[2],rawbuf[3],rawbuf[4]);
	
	// decode data in first word

	for (j = 0; j < 8; j++){
		preamble = preamble | (bitword[7 - j] << j);
	}
	for (j = 0; j < 6; j++){
		message_type = message_type + (bitword[13 - j] << j);
	}
	for (j = 0; j < 10; j++){
		station_id = station_id | (bitword[23 - j] << j);
	}

	loadword(rawbuf[5],rawbuf[6],rawbuf[7],rawbuf[8],rawbuf[9]);
	
	// decode data in second word

	for (j = 0; j < 13; j++){
		mod_zcount = mod_zcount | (bitword[12 - j] << j);
	}
	for (j = 0; j < 3; j++){
		seq_no = seq_no | (bitword[15 - j] << j);
	}
	for (j = 0; j < 5; j++){
		length = length | (bitword[20 - j] << j);
	}
	for (j = 0; j < 3; j++){
		health = health | (bitword[23 - j] << j);
	}

	// This is defined in the RTCM standard

	if (message_type == 0){
		message_type = 64;
	}
	
	// Added to help in debugging the output
	// Often some messages are lost in transmission

	fprintf(seqno,"%d\n",seq_no);
	
	// Set the correct filepointer so output is directed
	// to the correct place. Default is screen to flag an
	// error or warning
	
	switch(message_type){
		case 1: header = mtype1; break;
		case 2: header = mtype2; break;
		case 3: header = mtype3; break;
		case 6: header = mtype6; break;
		case 16: header = mtype16; break;
		case 18: header = mtype18; break;
		case 19: header = mtype19; break;
		default: header = stdout;
	}
	
	// Write header information in files
	
	if (bool_pre){
		fprintf(header,"%d ",preamble);
	}

	if (bool_mtype){
		fprintf(header,"%2.0f ",message_type);
	}
	
	if (bool_refid){
		fprintf(header,"%4.0f ",station_id);
	}

	if (bool_zcount){
		fprintf(header,"%6.1f ",mod_zcount*0.6);
	}

	if (bool_seqno){
		fprintf(header,"%d ",seq_no);
	}

	if (bool_length){
		fprintf(header,"%2.0f ",length);
	}

	if (bool_health){
		fprintf(header,"%d ",health);
	}
	
	// Call the correct function to decode the rest of the message.

	switch(message_type){
		case 1:	readmessage1(length); break;
		case 2: readmessage2(length); break;
		case 3: readmessage3(length); break;
		case 6: break;
		case 16: readmessage16(length); break;
		case 18: readmessage18(length); break;
		case 19: readmessage19(length); break;
		default: printf("WARNING: (l.%d) Unsupported message type detected! (type %d)\n",l,message_type);
	}
	
}

/*****************************************************************************
 * MAIN                                                                      *
 *                                                                           *
 * MAIN reads the command line arguments to determine which part of the data *
 * the user wants as output. If no arguments are given the help screen is    *
 * displayed. The input and output files are opened/created and the input    *
 * file is read line by line by searching for the CrLf symbol. If the size   *
 * of the message is not a multiple of five an error is flagged otherwise    *
 * the message is decode using the READHEADER function                       *
 *****************************************************************************/
 
void main(int argc,char *argv[]){
	
	// decode input arguments
	
	for (count = 0; count < argc; count++){
		if (strcmp(argv[count],"-cslip") == 0){
			bool_cslip = 1;
		}
		if (strcmp(argv[count],"-dprc") == 0){
			bool_dprc = 1;
		}
		if (strcmp(argv[count],"-drrc") == 0){
			bool_drrc = 1;
		}
		if (strcmp(argv[count],"-factor") == 0){
			bool_factor = 1;
		}
		if (strcmp(argv[count],"-freq") == 0){
			bool_freq = 1;
		}
		if (strcmp(argv[count],"-extime") == 0){
			bool_gnsstime = 1;
		}
		if (strcmp(argv[count],"-gpsglo") == 0){
			bool_gpsglo = 1;
		}
		if (strcmp(argv[count],"-health") == 0){
			bool_health = 1;
		}
		if (strcmp(argv[count],"-iod") == 0){
			bool_iod = 1;
		}
		if (strcmp(argv[count],"-length") == 0){
			bool_length = 1;
		}
		if (strcmp(argv[count],"-mpe") == 0){
			bool_mpe = 1;
		}
		if (strcmp(argv[count],"-msg") == 0){
			bool_msg = 1;
		}
		if (strcmp(argv[count],"-mtype") == 0){
			bool_mtype = 1;
		}
		if (strcmp(argv[count],"-multi") == 0){
			bool_multi = 1;
		}
		if (strcmp(argv[count],"-pca") == 0){
			bool_pcind = 1;
		}
		if (strcmp(argv[count],"-phase") == 0){
			bool_phase = 1;
		}
		if (strcmp(argv[count],"-pr") == 0){
			bool_pr = 1;
		}
		if (strcmp(argv[count],"-prc") == 0){
			bool_prc = 1;
		}
		if (strcmp(argv[count],"-pre") == 0){
			bool_pre = 1;
		}
		if (strcmp(argv[count],"-prn") == 0){
			bool_prn = 1;
		}
		if (strcmp(argv[count],"-dqual") == 0){
			bool_qual = 1;
		}
		if (strcmp(argv[count],"-refid") == 0){
			bool_refid = 1;
		}
		if (strcmp(argv[count],"-rrc") == 0){
			bool_rrc = 1;
		}
		if (strcmp(argv[count],"-seqno") == 0){
			bool_seqno = 1;
		}
		if (strcmp(argv[count],"-smooth") == 0){
			bool_smooth = 1;
		}
		if (strcmp(argv[count],"-udre") == 0){
			bool_udre = 1;
		}
		if (strcmp(argv[count],"-xyz") == 0){
			bool_xyz = 1;
		}
		if (strcmp(argv[count],"-zcount") == 0){
			bool_zcount = 1;
		}
		if (strcmp(argv[count],"-p00") == 0){
			d29star = 0;
			d30star = 0;
		}
		if (strcmp(argv[count],"-p01") == 0){
			d29star = 0;
			d30star = 1;
		}
		if (strcmp(argv[count],"-p10") == 0){
			d29star = 1;
			d30star = 0;
		}
		if (strcmp(argv[count],"-p11") == 0){
			d29star = 1;
			d30star = 1;
		}
		if (strcmp(argv[count],"-default") == 0){
			bool_dprc = 1;
			bool_drrc = 1;
			bool_prc = 1;
			bool_prn = 1;
			bool_rrc = 1;
			bool_xyz = 1;
			bool_zcount = 1;
			bool_msg = 1;
			bool_phase = 1;
			bool_pr = 1;
		}
	}
	
	// if no arguments are given then display help screen and terminate execution
	
	if  (!(bool_cslip || bool_dprc || bool_drrc || bool_factor || bool_freq || bool_gnsstime || bool_gpsglo || bool_health || bool_iod || bool_length || bool_mpe || bool_msg || bool_mtype || bool_multi || bool_pcind || bool_phase || bool_pr || bool_prc || bool_pre || bool_prn || bool_qual || bool_refid || bool_rrc || bool_seqno || bool_smooth || bool_udre || bool_xyz || bool_zcount)){
		printf("\nUSAGE: rtcm2asc.exe [-flags [-flags... ]]\n\n");
		
		printf("This program converts the binary RTCM messages to ASCII format.\n");
		printf("The program requires the RTCM file to be named \"rtcmout.bin\" and placed in the\n");
		printf("same folder as the .exe file\n\n");
		
		printf("ASCII text files named \"mtype1.txt\",\"mtype2.txt\" etc. are then generated by\n");
		printf("the program and can be easily loaded into eg. Matlab\n");
		printf("Currently only message types: 1,2,3,6,16,18 and 19 are supported\n\n");
		printf("The amount of data stored in the ASCII files are controlled by using\n");
		printf("command-line arguments:\n\n");
		
		printf("\"-cslip\":   Cumulative loss of continuity indicator [0 - 31]\n");
		printf("\"-dprc\":    Delta Pseudorange corrections [+-655.34 or +-10485.44 m]\n");
		printf("\"-drrc\":    Delta Range-Rate Correction [+-0.254 or +- 4.064 m/s]\n");
		printf("\"-dqual\":   Data quality indicator\n");
		printf("\"-extime\":  Extended GNSS time [0-599999 us]\n");
		printf("\"-factor\":  Scale factor bit [0|1]\n");
		printf("\"-freq\":    Frequency indicator [0-3] 0 = L1, 2 = L2\n");
		printf("\"-gpsglo\":  GPS/GLONASS indicator [0|1] 0 = GPS, 1 = GLONASS\n\n");

		printf("Press any key to continue...");

		getch();
		
		printf("\n\n\"-health\":  Health of the reference station [0-7]\n");
		printf("\"-iod\":     Issue-of-data of the satellites [0-256]\n");
		printf("\"-length\":  Number of words in message containing data\n");
		printf("\"-mpe\":     Multipath error [0-15]\n");
		printf("\"-msg\":     Message from reference station, max 90 characters long\n");
		printf("\"-mtype\":   Message type [1-64]\n");
		printf("\"-multi\":   Multi-message indicator [0|1]\n");
		printf("\"-pca\":     P-code or C/A-code indicator [0|1] 0 = C/A, 1 = P\n");
		printf("\"-phase\":   Uncorrected carrier phases [+-8,388,608 cycles]\n");
		printf("\"-pr\":      Uncorrected pseudoranges [0-85,899,345.90 m]\n");
		printf("\"-prc\":     Pseudorange corrections [+-655.34 or +-10,485.44 m]\n");
		printf("\"-pre\":     Preamble [102]\n");
		printf("\"-prn\":     Satellite ID [1-32]\n");
		printf("\"-pXY\":     Enter initial parity values for bits 29 and 30. X [0|1], Y [0|1]\n");
		printf("\"-refid\":   Reference station ID [0-1023]\n");
		printf("\"-rrc\":     Range-Rate Correction [+-0.254 or +-4.064 m/s]\n");
		printf("\"-seqno\":   Sequence number of the message [0-7]\n");
		printf("\"-smooth\":  Smoothing interval [0-3]\n");
		printf("\"-udre\":    User Differential Range Error [0-3]\n");
		printf("\"-xyz\":     ECEF position (XYZ) of the ref. station [0-21,474,836.47 m]\n");
		printf("\"-zcount\":  Modified z-count [0-3599.4 sec]\n\n");

		printf("\"-default\": Uses the most common of the above parameters\n");
		exit(0);
	}
	
	// else open input files for reading and writing
	
	if ((rtcmbin  = fopen("rtcmout.bin","rb")) == NULL){
      printf("ERROR opening input file\n" );
	  exit(1);
	}

	if ((mtype1  = fopen("mtype1.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype1.txt'\n" );
	  exit(1);
	}

	if ((mtype2  = fopen("mtype2.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype2.txt'\n" );
	  exit(1);
	}

	if ((mtype3  = fopen("mtype3.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype3.txt'\n" );
	  exit(1);
	}

	if ((mtype6  = fopen("mtype6.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype6.txt'\n" );
	  exit(1);
	}

	if ((mtype16  = fopen("mtype16.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype16.txt'\n" );
	  exit(1);
	}
	
	if ((mtype18  = fopen("mtype18.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype18.txt'\n" );
	  exit(1);
	}

	if ((mtype19  = fopen("mtype19.txt","wt")) == NULL){
      printf("ERROR opening output file 'mtype19.txt'\n" );
	  exit(1);
	}
	
	if ((seqno  = fopen("seqno.txt","wt")) == NULL){
      printf("ERROR opening output file 'seqno.txt'\n" );
	  exit(1);
	}

	// for each line in the input use readheader
	// to decode the contents

	l = 0;
	while (!feof(rtcmbin)){
		k = 0;
		d = 0;
		rawbuf[0] = '\0';
		while ((d == 0) && (k < 255)){
			ch = fgetc(rtcmbin);
			
			// other termination strings could be used here
			
			if (ch == 0x0D){
				if (0x0A == fgetc(rtcmbin)) {
					d = 1;
					l++;
				}
			}
			else {
				
				// if not termination string then read into buffer
				rawbuf[k] = ch;
				k++;
			}
		}
		rawbuf[k] = '\0';

		div_result = div(k,5);

		// check if message size is a multiple of five

		if (div_result.rem != 0){
			printf("ERROR: (l.%d) Input does not follow word boundaries!\n\n",l);
		}

		readheader();
	}

	// finish by closing all open files

	_fcloseall();
}