#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define MAX_NO_RE   100
#define MAX_NO_AA   64
#define FILE_NAME_SIZE 12
#define MAX_INPUT_LEN  256
#define AA_PER_LINE  35

typedef struct {
	char f_aa[20];
	char s_aa[20];
	char t_aa[20];
} RF;

typedef struct {
	char name[45];
	char na[7];
} RE;

typedef struct {
	char nucleic_acid[4];
	char aa;
} AA;

typedef struct {
	int pos;
	int number;
	char name[45];
}  OUTPUT;

RF f_rf[MAX_NO_RE], s_rf[MAX_NO_RE], t_rf[MAX_NO_RE];
RE res_enzyme[MAX_NO_RE];
AA amino_acid[MAX_NO_AA];
OUTPUT out[200];
char base[4] = { 'A', 'C', 'G', 'T' };

int nre, naa, nout;

int IsChIn(char *str, char c);
void main(int argc, char *argv[]);
void DisplayReTable(FILE *fp);

/******************************************************************************
*                                                                             *
*   ReadDataBase_RE:    Reads the amino acid sequence motifs compatible with  *
*                       the Restriction Enzyme recognition sequences from a   *
*                       file.                                                 *
*                                                                             *
*   Input:              fname. File containing the amino acid sequence motifs *
*                       for all Restriction Enzymes.                          *
*                                                                             *
*   Output:		None.                                                     *
*                                                                             *
*   Notes:              This function reads the nucleic acid sequence  for    *
*                       for all restriction enzymes and generates the amino   *
*               	    acid motifs from each of the reading frame for each of*
*                       the restriction enzymes. This information is used later
*                       to determine potential sites for mutation in a given  *
*                       sequence of amino acids.                              *
*                                                                             *
******************************************************************************/

ReadDataBase_RE(char *fname)
{
	int i, j;
	FILE *fp;
	char c, re[7];
	char name[50];

	/*  Open the file for reading the amino acid sequence of restriction enzymes */

	if ((fp = fopen(fname, "r")) == (FILE *)NULL)
	{
		printf("Error opening DataBase file %s\n", fname);
		exit(-1);
	}

	/* Initialize the reading frames */

	for (i = 0; i < MAX_NO_RE; i++)
	{
		for (j = 1; j < 20; j++)
		{
			f_rf[i].f_aa[j] = '\0';
			f_rf[i].s_aa[j] = '\0';
			f_rf[i].t_aa[j] = '\0';
			s_rf[i].f_aa[j] = '\0';
			s_rf[i].s_aa[j] = '\0';
			s_rf[i].t_aa[j] = '\0';
			t_rf[i].f_aa[j] = '\0';
			t_rf[i].s_aa[j] = '\0';
			t_rf[i].t_aa[j] = '\0';
		}
	}

	nre = 0;

	/* read the amino acid sequence and name of each restriction enzyme */
	c = fgetc(fp);
	while (!feof(fp))
	{
		i = 0;
		while (c != ' ')
		{
			re[i++] = c;
			c = fgetc(fp);
		}
		re[i] = '\0';
		for (; c == ' ';)
			c = fgetc(fp);

		i = 0;
		while ((c != '\n') && (c != -1))
		{
			name[i++] = c;
			c = fgetc(fp);
		}
		name[i] = '\0';

		for (; c == '\n';)
			c = fgetc(fp);

		/* Generate the First reading frame */
		First_Rf(re, nre);

		/* Generate the Second reading frame */
		Second_Rf(re, nre);

		/* Generate the Third reading frame */
		Third_Rf(re, nre);

		/* Store the name and amino acid sequence of the RE */
		strcpy(res_enzyme[nre].name, name);
		strcpy(res_enzyme[nre].na, re);

		nre++;
	}
	fclose(fp);
}


/******************************************************************************
*                                                                             *
*   First_Rf:	Generates the first reading frame for a Restriction Enzyme    *
*                                                                             *
*   Input:	nucliec acid sequence of a RE.                                *
*                                                                             *
*   Output:		None.                                                 *
*                                                                             *
*   Notes:	The first reading frame contains two amino acids. A single    *
*		letter code corresponding to these amino acids are stored in  *
*		the first reading frame.                                      *
*                                                                             *
******************************************************************************/

First_Rf(char *str, int n)
{
	int i;
	char str1[4], str2[4];

	strncpy(str1, str, 3);
	str1[3] = '\0';
	strncpy(str2, &str[3], 3);
	str2[3] = '\0';

	/* Determine the first amino acid */

	for (i = 0; i < MAX_NO_AA; i++)
	{
		if ((strcmp(amino_acid[i].nucleic_acid, str1)) == 0)
		{
			f_rf[n].f_aa[0] = amino_acid[i].aa;
			f_rf[n].f_aa[1] = '\0';
			break;
		}
	}

	/* Determine the second amino acid */

	for (i = 0; i < MAX_NO_AA; i++)
	{
		if ((strcmp(amino_acid[i].nucleic_acid, str2)) == 0)
		{
			f_rf[n].s_aa[0] = amino_acid[i].aa;
			f_rf[n].s_aa[1] = '\0';
			break;
		}
	}
}

/******************************************************************************
*                                                                             *
*   Second_Rf:	Generates the second reading frame for a Restriction Enzyme   *
*                                                                             *
*   Input:	nucliec acid sequence of a RE.                                *
*                                                                             *
*   Output:	None.                                                         *
*                                                                             *
*   Notes:	The second first reading frame contains three amino acids. A  *
*		single letter code corresponding to these amino acids are     *
*		stored in the second reading frame.                           *
*                                                                             *
******************************************************************************/

Second_Rf(char *str, int n)
{
	int i, j;
	char str1[4];

	/* Determine the first amino acid */

	for (i = 0, j = 0; i < MAX_NO_AA; i++)
	{
		if ((str[0] == amino_acid[i].nucleic_acid[1]) &&
			(str[1] == amino_acid[i].nucleic_acid[2]))
		{
			if (!IsChIn(s_rf[n].f_aa, amino_acid[i].aa))
				s_rf[n].f_aa[j++] = amino_acid[i].aa;
		}
	}
	strncpy(str1, &str[2], 3);
	str1[3] = '\0';

	/* Determine the second amino acid */

	for (i = 0; i < MAX_NO_AA; i++)
	{
		if ((strcmp(str1, amino_acid[i].nucleic_acid)) == 0)
		{
			s_rf[n].s_aa[0] = amino_acid[i].aa;
			break;
		}
		s_rf[n].s_aa[1] = '\0';
	}

	/* Determine the third amino acid */

	for (i = 0, j = 0; i < MAX_NO_AA; i++)
	{
		if (str[5] == amino_acid[i].nucleic_acid[0])
		{
			if (!IsChIn(s_rf[n].t_aa, amino_acid[i].aa))
				s_rf[n].t_aa[j++] = amino_acid[i].aa;
		}
	}
}

/******************************************************************************
*                                                                             *
*   Third_Rf:	Generates the third  reading frame for a Restriction Enzyme   *
*                                                                             *
*   Input:	nucliec acid sequence of a RE.                                *
*                                                                             *
*   Output:	None.                                                         *
*                                                                             *
*   Notes:	The third  first reading frame contains three amino acids. A  *
*		single letter code corresponding to these amino acids are     *
*		stored in the third  reading frame.                           *
*                                                                             *
******************************************************************************/


Third_Rf(char *str, int n)
{
	int i, j;
	char str1[4];

	/* Determine the second amino acid */

	for (i = 0, j = 0; i < MAX_NO_AA; i++)
	{
		if ((str[4] == amino_acid[i].nucleic_acid[0]) &&
			(str[5] == amino_acid[i].nucleic_acid[1]))
		{
			if (!IsChIn(t_rf[n].t_aa, amino_acid[i].aa))
				t_rf[n].t_aa[j++] = amino_acid[i].aa;
		}
	}

	strncpy(str1, &str[1], 3);
	str1[3] = '\0';

	/* Determine the first  amino acid */

	for (i = 0; i < MAX_NO_AA; i++)
	{
		if ((strcmp(str1, amino_acid[i].nucleic_acid)) == 0)
		{
			t_rf[n].s_aa[0] = amino_acid[i].aa;
			break;
		}
		t_rf[n].s_aa[1] = '\0';
	}

	/* Determine the third  amino acid */

	for (i = 0, j = 0; i < MAX_NO_AA; i++)
	{
		if (str[0] == amino_acid[i].nucleic_acid[2])
		{
			if (!IsChIn(t_rf[n].f_aa, amino_acid[i].aa))
				t_rf[n].f_aa[j++] = amino_acid[i].aa;
		}
	}
}

/******************************************************************************
*                                                                             *
*   ReadDataBase_AA: 	Reads the nucliec acid sequence for all the amino     *
*			acids from a file.                                    *
*                                                                             *
*   Input:		fname. File containing the nucliec acid sequence for  *
*			all amino acids.                                      *
*                                                                             *
*   Output:		None.                                                 *
*                                                                             *
*   Notes:		This function reads the nucliec acid sequence and the *
*			one letter code for all the amino acids.              *
*                                                                             *
******************************************************************************/

ReadDataBase_AA(char *fname)
{
	int i;
	FILE *fp;
	char c;

	if ((fp = fopen(fname, "r")) == (FILE *)NULL)
	{
		printf("Error opening DataBase file %s\n", fname);
		exit(-1);
	}

	naa = 0;

	c = fgetc(fp);
	while (!feof(fp))
	{

		i = 0;
		while (c != ' ')
		{
			amino_acid[naa].nucleic_acid[i++] = c;
			c = fgetc(fp);
		}
		for (; c == ' ';)
			c = fgetc(fp);

		amino_acid[naa].aa = c;
		naa++;
		c = fgetc(fp);
		for (; c == '\n';)
			c = fgetc(fp);
	}
	fclose(fp);
}

/******************************************************************************
*                                                                             *
*   DisplayReTable: 	Writes the reading frame corresponding to all the     *
*			Restriction Enzymes onto a file.                      *
*                                                                             *
*   Input:		fp. file pointer for writing.                         *
*                                                                             *
*   Output:		None.                                                 *
*                                                                             *
*   Notes:		This function can be used to generate the table       *
*			containg the reading frames for Restriction enzymes.  *
*                                                                             *
******************************************************************************/

void DisplayReTable(FILE *fp)
{
	int i, j;


	for (i = 0; i < nre; i++)
	{
		fprintf(fp, "%s", res_enzyme[i].name);
		for (j = strlen(res_enzyme[i].name); j < 45; j++)
			fprintf(fp, " ");
		fprintf(fp, "%s", res_enzyme[i].na);
		fprintf(fp, "  ");
		fprintf(fp, "%c ", f_rf[i].f_aa[0]);
		fprintf(fp, "%c  ", f_rf[i].s_aa[0]);

		fprintf(fp, "%s", s_rf[i].f_aa);
		for (j = strlen(s_rf[i].f_aa); j < 5; j++)
			fprintf(fp, " ");
		fprintf(fp, " %c ", s_rf[i].s_aa[0]);
		fprintf(fp, "%s", s_rf[i].t_aa);
		for (j = strlen(s_rf[i].t_aa); j < 8; j++)
			fprintf(fp, " ");


		fprintf(fp, "%s", t_rf[i].f_aa);
		for (j = strlen(t_rf[i].f_aa); j < 15; j++)
			fprintf(fp, " ");
		fprintf(fp, " %c ", t_rf[i].s_aa[0]);
		fprintf(fp, "%s", t_rf[i].t_aa);

		fprintf(fp, "\n");
	}

}

int IsChIn(char *str, char c)
{
	int i;

	for (i = 0; str[i]; i++)
		if (str[i] == c)
			return(1);

	return(0);
}

void main(int argc, char *argv[])
{
	char aa_database[FILE_NAME_SIZE];
	char re_database[FILE_NAME_SIZE];
	char fname[FILE_NAME_SIZE];
	FILE *fp;


	printf("Program to generate Restriction Enzyme Table\n");

	strcpy(aa_database, "dbase1");
	ReadDataBase_AA(aa_database);

	strcpy(re_database, "dbase2");
	ReadDataBase_RE(re_database);

	if (argc > 1)
	{
		if ((fp = fopen(argv[1], "w+")) == (FILE *)NULL)
			fp = stdout;
	}
	else
		fp = stdout;

	DisplayReTable(fp);

}
