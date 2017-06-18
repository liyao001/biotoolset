
/****************************************************************************
*                                                                           *
*       Program (SILMUT) for recognizing the potential mutation sites       *
*       in an amino acid sequence.                                          *
*       Author: K Vijayananda, Department of Computer Science,              *
*       University of Maryland, College Park, MD 20740                      *
****************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>


#define MAX_NO_RE   100
#define MAX_NO_AA   64
#define MAX_NO_NA   4
#define MAX_MS   200
#define FILE_NAME_SIZE 12
#define MAX_INPUT_LEN  256
#define AA_PER_LINE  35
#define LINELEN 80
typedef struct
{
    char f_aa[20];
    char s_aa[20];
    char t_aa[20];
} RF;

typedef struct
{
    char name[45];
    char na[7];
} RE;

typedef struct
{
    char nucleic_acid[4];
    char aa;
} AA;

typedef struct
{
    int pos;
    int number;
    char name[45];
    char na[7];
}  OUTPUT;

RF f_rf[MAX_NO_RE], s_rf[MAX_NO_RE], t_rf[MAX_NO_RE];
RE res_enzyme[MAX_NO_RE];
AA amino_acid[MAX_NO_AA];
OUTPUT out[MAX_MS];
char base[5] = { 'A', 'C', 'G', 'T', '\0' };
char valid_aa[MAX_NO_AA + 1];
int nre, naa, nout;

int IsChIn(char *str, char c);
int ConvertNAToAA(char *in_str, char **aa, int option);
int Duplicate(char *str[], int n);
void  main(int argc, char *argv[]);

/******************************************************************************
*                                                                             *
*   ReadDataBase_RE:    Reads the nucleic acid recognition sequence for       *
*                       the Restriction Enzymes from a file.                  *
*                                                                             *
*                                                                             *
*   Input:              fname. File containing the nucleic acid sequences     *
*                       for the Restriction Enzymes.                          *
*                                                                             *
*   Output:             None.                                                 *
*                                                                             *
*   Notes:              This function reads the nucleic acid sequences for    *
*                       for the restriction enzymes and generates the amino   *
*                       acid motifs from each of the reading frame for each of*
*                       the restriction enzymes. This information is used     *
*                       later to determine potential sites for mutation in a  *
*                       given or derived sequence of amino acids.             *
*                                                                             *
******************************************************************************/

ReadDataBase_RE(char *fname)
{
    int i, j;
    FILE *fp;
    char c, re[7];
    char name[50];

    /*  Open the file for reading the nucleic acid sequence of restriction enzymes */
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

    /* read the nucleic acid sequence and name of each restriction enzyme */
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
*   First_Rf:   Generates the amino acid motifs from the first reading frame  *
*               for a Restriction Enzyme                                      *
*                                                                             *
*   Input:      nucleic acid sequence of a RE.                                *
*                                                                             *
*   Output:     None.                                                         *
*                                                                             *
*   Notes:      The first reading frame contains two amino acids for a        *
*               six-base recognition sequence.  A single letter code          *
*               corresponding to these amino acids are stored in              *
*               the first reading frame.                                      *
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
*   Second_Rf: Generates the amino acid motifs from the second reading frame  *
*              for a Restriction Enzyme                                       *
*   Input:     nucleic acid sequence of a RE.                                 *
*                                                                             *
*   Output:    None.                                                          *
*                                                                             *
*   Notes:     The second reading frame contains three amino acids for a      *
*              six base recognition sequence. A single letter code            *
*              corresponding to these amino acids are stored.                 *
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
*   Third_Rf: Generates the amino acid motifs from the third reading frame    *
*             for a Restriction Enzyme                                        *
*   Input:    nucleic acid sequence of a RE.                                  *
*                                                                             *
*   Output:   None.                                                           *
*                                                                             *
*   Notes:    The third reading frame contains three amino acids for a        *
*             six base recognition sequence. A single letter code             *
*             corresponding to these amino acids are stored.                  *
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
*   ReadDataBase_AA:    Reads the nucleic acid codons for all the amino       *
*                       acids from a file.                                    *
*                                                                             *
*   Input:              fname. File containing the nucleic acid codons for    *
*                       all amino acids.                                      *
*                                                                             *
*   Output:             None.                                                 *
*                                                                             *
*   Notes:              This function reads the nucleic acid codons and the   *
*                       corresponding one letter code  all the amino acids.   *
*                                                                             *
******************************************************************************/

ReadDataBase_AA(char *fname)
{
    int i, j;
    FILE *fp;
    char c;

    /*  Open the file for reading the nucleic acid codons and the corresponding
    one letter code  for all amino acids */

    if ((fp = fopen(fname, "r")) == (FILE *)NULL)
    {
        printf("Error opening DataBase file %s\n", fname);
        exit(-1);
    }

    naa = 0;
    j = 0;
    valid_aa[0] = '\0';

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
        if (!IsChIn(valid_aa, c))
            valid_aa[j++] = c;

        naa++;
        c = fgetc(fp);
        for (; c == '\n';)
            c = fgetc(fp);
    }
    valid_aa[j] = '\0';
    fclose(fp);
}

/*******************************************************************************
*                                                                              *
*   DisplayReTable:   Writes the amino acid motifs from each of the reading    *
*                     frame corresponding to all the Restriction Enzymes       *
*                     onto a file.                                             *
*                                                                              *
*   Input:            fp. file pointer for writing.                            *
*                                                                              *
*   Output:           None.                                                    *
*                                                                              *
*   Notes:            This function generates the table containing the amino   *
*                     acid motifs for each of the reading frames corresponding *
*                     to all the restriction enzymes.                          *
*                                                                              *
*******************************************************************************/

DisplayReTable(FILE *fp)
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

/******************************************************************************
*                                                                             *
*   ConverNAToAA:   converts a nucleic acid sequence into its corresponding   *
*                   amino acid sequence.                                      *
*                                                                             *
*   Input: in_str:  nucleic acid sequence.                                    *
*                   aa:  pointer to string for storing the amino acid         *
*                   sequence.                                                 *
*                   option:  method for conversion.                           *
*                                                                             *
*   Output:         returns the number of amino acid sequences generated.     *
*                                                                             *
*   Notes:          This function converts a nucleic acid sequence into       *
*                   its corresponding amino acid sequence. If the input       *
*                   sequence is not a multiple of 3, the option parameter     *
*                   determines the recovery action. One recovery action       *
*                   is to determine all possible amino acid sequence for      *
*                   the given input. If the input sequence is one short,      *
*                   4 nucleic acid sequences are generated. If it is 2 short, *
*                   then 16 nucleic acid sequence are generated.              *
*                                                                             *
******************************************************************************/

int ConvertNAToAA(char *in_str, char **aa, int option)
{
    int i, j, k, n, m, len, count;
    char temp[4], *str;

    len = strlen(in_str);
    str = (char *)calloc(len + 4, sizeof(char));
    strcpy(str, in_str);

    switch (option)
    {
    case 0:

        aa[0] = (char *)calloc(len + 1, sizeof(char));

        for (i = 0, k = 0; i < len; i += 3)
        {
            temp[0] = str[i];
            temp[1] = str[i + 1];
            temp[2] = str[i + 2];
            temp[3] = '\0';

            for (j = 0; j < naa; j++)
            {
                if (!strcmp(amino_acid[j].nucleic_acid, temp))
                {
                    aa[0][k++] = amino_acid[j].aa;
                    break;
                }
            }
        }
        aa[0][k] = '\0';
        count = 1;
        break;

    case 2:

        str[len + 1] = '\0';
        for (n = 0; n < 4; n++)
        {
            str[len] = base[n];
            aa[n] = (char *)calloc(len + 4, sizeof(char));

            for (i = 0, k = 0; i < len + 1; i += 3)
            {
                temp[0] = str[i];
                temp[1] = str[i + 1];
                temp[2] = str[i + 2];
                temp[3] = '\0';

                for (j = 0; j < naa; j++)
                {
                    if (!strcmp(amino_acid[j].nucleic_acid, temp))
                    {
                        aa[n][k++] = amino_acid[j].aa;
                        break;
                    }
                }
            }
            aa[n][k] = '\0';
        }
        count = 4;
        break;


    case 1:

        str[len + 2] = '\0';
        for (n = 0; n < 4; n++)
        {
            str[len] = base[n];
            for (m = 0; m < 4; m++)
            {
                str[len + 1] = base[m];
                aa[n * 4 + m] = (char *)calloc(len + 4, sizeof(char));

                for (i = 0, k = 0; i < len + 2; i += 3)
                {
                    temp[0] = str[i];
                    temp[1] = str[i + 1];
                    temp[2] = str[i + 2];
                    temp[3] = '\0';

                    for (j = 0; j < naa; j++)
                    {
                        if (!strcmp(amino_acid[j].nucleic_acid, temp))
                        {
                            aa[n * 4 + m][k++] = amino_acid[j].aa;
                            break;
                        }
                    }
                }
                aa[n * 4 + m][k] = '\0';
            }
        }
        count = 16;
        break;

    default:
        count = 0;

    }
    return(count);

}

/*******************************************************************************
*                                                                              *
*   ScanForRE:  For a given string of amino acids, it recognizes the potential *
*               sites for mutation and the restriction enzymes that can be     *
*               introduced at this particular site.                            *
*                                                                              *
*   Input:      string of amino acids.                                         *
*                                                                              *
*   Output:     None.                                                          *
*                                                                              *
*   Notes:      This function scans the amino acid sequence and checks for the *
*               existence of amino acid motifs obtained from each of the       *
*               reading frames (first, second and third) of all the            *
*               Restriction enzymes in the input amino acid sequence.          *
*               The name of the restriction enzyme and the location of the     *
*               site in input sequence are stored in anarray.                  *
*                                                                              *
*******************************************************************************/
ScanForRE(char *str)
{

    int i, j, k, len;

    len = strlen(str);

    for (i = 0, k = 0; i < len - 2; i++)
    {
        /* Check for first  reading frame */
        for (j = 0; j < nre; j++)
        {
            if ((IsChIn(f_rf[j].f_aa, str[i])) && (IsChIn(f_rf[j].s_aa, str[i + 1])))
            {
                out[k].pos = i;
                out[k].number = 2;
                strcpy(out[k].name, res_enzyme[j].name);
                strcpy(out[k].na, res_enzyme[j].na);
                k++;
            }
        }

        /* Check for second reading frame */
        for (j = 0; j < nre; j++)
        {
            if ((IsChIn(s_rf[j].f_aa, str[i])) && (IsChIn(s_rf[j].s_aa, str[i + 1]))
                    && (IsChIn(s_rf[j].t_aa, str[i + 2])))
            {
                out[k].pos = i;
                out[k].number = 3;
                strcpy(out[k].name, res_enzyme[j].name);
                strcpy(out[k].na, res_enzyme[j].na);
                k++;
            }
        }
        /* Check for third  reading frame */

        for (j = 0; j < nre; j++)
        {
            if ((IsChIn(t_rf[j].f_aa, str[i])) && (IsChIn(t_rf[j].s_aa, str[i + 1]))
                    && (IsChIn(t_rf[j].t_aa, str[i + 2])))
            {
                out[k].pos = i;
                out[k].number = 3;
                strcpy(out[k].name, res_enzyme[j].name);
                strcpy(out[k].na, res_enzyme[j].na);
                k++;
            }
        }

    }
    for (j = 0; j < nre; j++)
    {
        if ((IsChIn(f_rf[j].f_aa, str[i])) && (IsChIn(f_rf[j].s_aa, str[i + 1])))
        {
            out[k].pos = i;
            out[k].number = 2;
            strcpy(out[k].name, res_enzyme[j].name);
            strcpy(out[k].na, res_enzyme[j].na);
            k++;
        }
    }
    nout = k;
}

/*******************************************************************************
*                                                                              *
*   GetNum:     Get a number form a given file.                                *
*                                                                              *
*   Input:      File Pointer                                                   *
*                                                                              *
*   Output:     -1 if EOF is reached                                           *
*               0  if invalid input                                            *
*               n  the number that is input                                    *
*   Notes:      All trailing blanks are ignored.                               *
*               ignored.
*                                                                              *
*******************************************************************************/

int GetNum(FILE *fp)
{
    int i, j, c, opt;
    char str[20];

    c = fgetc(fp);
    i = 0;
    while ((c != '\n') && (c != EOF))
    {
        str[i++] = c;
        c = fgetc(fp);
    }
    str[i] = '\0';

    if (c == EOF)
        return(EOF);

    if (i == 0)
        return(0);

    j = 0;
    opt = 0;
    while (j < i)
    {
        if (isdigit(str[j]))
        {
            opt = opt * 10 + (str[j] - '0');
            j++;
        }
        else if (str[j] == ' ')
        {
            while (str[j] == ' ')
                j++;
        }
        else
            return(0);
    }

    return(opt);

}

int Check_Input(char *str, int opt)
{
    int i, j, len;
    char c;

    len = strlen(str);
    if (len == 0)
        return(0);

    if (str[len - 1] == '\n')
        str[len - 1] = '\0';


    len = strlen(str);
    for (i = 0, j = 0; i < len; i++)
    {
        if (str[i] != ' ')
            str[j++] = str[i];
    }
    str[j] = '\0';

    for (i = 0; str[i]; i++)
    {
        c = toupper(str[i]);
        str[i] = c;
    }

    if (opt == 1)
    {
        for (i = 0; str[i]; i++)
        {
            if (!IsChIn(valid_aa, str[i]))
                return(0);
        }
    }
    else if (opt == 2)
    {
        for (i = 0; str[i]; i++)
        {
            if (!IsChIn(base, str[i]))
                return(0);
        }
    }

    return(1);
}

/******************************************************************************
*                                                                             *
*   PrintResult:        print the result indicating the mutation site and the *
*                       name of the restriction enzyme that can be introduced *
*                       at this site                                          *
*                                                                             *
*   Input:              string of amino acids and the file for output.        *
*                                                                             *
*   Output:             None.                                                 *
*                                                                             *
*   Notes:              This function prints out the amino acid sequence and  *
*                       the position in the string and the name of the        *
*                       restriction enzyme/s that match with the amino acid   *
*                       motif in the string.                                  *
*                                                                             *
*                                                                             *
******************************************************************************/
PrintResult(str, fp)
char *str;
FILE *fp;

{
    int i, j, k, l, m, len, pos;

    len = strlen(str);
    fprintf(fp, "\n\n-----------------------------------------------------------------------\n");

    if (nout <= 0)
    {
        for (i = 0; i < len; i++)
            fprintf(fp, "%c", str[i]);
        fprintf(fp, "\n");

        fprintf(fp, "No site in the input string can be replaced with Restriction Enzymes\n");
        return 0;
    }

    pos = 0;
    i = 0;

    while (i < nout)
    {
        j = i;
        while (j < nout)
        {
            len = out[j].pos - pos + strlen(out[j].name);
            if (len > LINELEN)
                break;
            else
                j++;
        }

        l = i;
        if (i == j)
        {
            len = out[i].pos;
            for (k = pos; k < len; k++)
                fprintf(fp, "%c", str[k]);

            pos += len;
        }
        else
        {
            if (j == nout)
                len = strlen(str) + 1;
            else
                len = out[j - 1].pos + out[j - 1].number;

            for (k = pos; k < len; k++)
                fprintf(fp, "%c", str[k]);
            fprintf(fp, "\n");

            k = pos;
            for (; i < j; i++)
            {
                for (; k < out[i].pos; k++)
                    fprintf(fp, " ");
                fprintf(fp, "%s", out[i].name);
                k += strlen(out[i].name);

                if ((i + 1 < j) && (k >= out[i + 1].pos))
                {
                    k = pos;
                    fprintf(fp, "\n");
                }
            }

            fprintf(fp, "\n\n");
            for (k = l; k < j; k++)
            {
                fprintf(fp, "Position in the input string: %d\n", out[k].pos - pos + 1);
                fprintf(fp, "Amino acid string at this position: ");
                if (out[k].number == 2)
                {
                    fprintf(fp, "%c", str[out[k].pos]);
                    fprintf(fp, "%c", str[out[k].pos + 1]);
                }
                else if (out[k].number == 3)
                {
                    fprintf(fp, "%c", str[out[k].pos]);
                    fprintf(fp, "%c", str[out[k].pos + 1]);
                    fprintf(fp, "%c", str[out[k].pos + 2]);
                }
                fprintf(fp, "\n");
                fprintf(fp, "Restriction Enzyme site/s that can be introduced at this position: ");
                fprintf(fp, "%s (%s)", out[k].name, out[k].na);
                fprintf(fp, "\n\n");
            }

            if (i < nout)
            {
                if ((out[i - 1].pos + out[i - 1].number) < out[i].pos)
                    pos += out[i - 1].pos + out[i - 1].number;
                else
                    pos += out[i - 1].pos;
            }

            if (j == nout)
                pos = strlen(str) - 1;
        }

        fprintf(fp, "\n");
    }

    len = strlen(str);
    if (pos < len - 1)
    {
        for (; pos < len; pos++)
            fprintf(fp, "%c", str[pos]);
    }
    fprintf(fp, "\n");

}

/******************************************************************************
*                                                                             *
*   Duplicate:      checks for a duplicate amino acid sequence.               *
*                                                                             *
*   Input:          string of amino acids and the file for output.            *
*                                                                             *
*   Output:   1:    if there is a duplicate                                   *
*             0:    otherwise                                                 *
*                                                                             *
*   Notes:          When the input nucleic acid sequence is not a multiple of *
*                   3, duplicate amino acid sequences might be generated.     *
*                   This function eliminates the analysis of duplicate        *
*                   amino acid sequences.                                     *
*                                                                             *
******************************************************************************/


int Duplicate(char *str[], int n)
{
    int i, res = 0;

    for (i = 0; i < n; i++)
    {
        if (!strcmp(str[i], str[n]))
            res = 1;
    }
    return(res);
}

void  main(argc, argv)
int argc;
char *argv[];

{
    char aa_database[FILE_NAME_SIZE];
    char re_database[FILE_NAME_SIZE];
    char input_str[MAX_INPUT_LEN], *aa_str[64];
    int option, i, c, len, n;
    FILE *res, *in;

    res = stdout;
    in = stdin;

    i = 1;
    while (i < argc)
    {
        if (!strcmp(argv[i], "-i"))
        {
            i++;
            if ((in = fopen(argv[i], "r")) == (FILE *)NULL)
                in = stdin;
        }
        else if (!strcmp(argv[i], "-o"))
        {
            i++;
            if ((res = fopen(argv[i], "w")) == (FILE *)NULL)
                res = stdin;
        }
        else
        {
            fprintf(stderr, "Invalid Option %s \n", argv[i]);
            fprintf(stderr, "Usage %s [-i <infile> -o <outfile>]\n", argv[0]);
            exit(-1);
        }
        i++;

    }

    strcpy(aa_database, "dbase1");
    ReadDataBase_AA(aa_database);

    strcpy(re_database, "dbase2");
    ReadDataBase_RE(re_database);

    while (1)
    {
        if (in == stdin)
        {
            printf("1:  Amino Acid Sequence.\n");
            printf("2:  Nucleic Acid Sequence.\n");
            printf("3:  Quit \n");
            printf("Enter number for the type of input sequence or 3 to quit: ");
        }

        option = GetNum(in);

        if (option == 3)
            break;

        if ((option == 1) || (option == 2))
        {
            if (in == stdin)
                printf("Enter the Input Sequence\n");

            i = 0;
            c = fgetc(in);
            while ((c != '\n') && (c != EOF))
            {
                input_str[i++] = c;
                c = fgetc(in);
            }
            input_str[i] = '\0';


            if (Check_Input(input_str, option))
            {
                len = strlen(input_str);

                if (option == 2)
                {
                    n = ConvertNAToAA(input_str, aa_str, (len % 3));
                    for (i = 0; i < n; i++)
                    {
                        /* Avoids analysis of duplicate amino acid sequences. */

                        if (!Duplicate(aa_str, i))
                        {
                            ScanForRE(aa_str[i]);
                            PrintResult(aa_str[i], res);
                        }
                    }
                }
                else
                {
                    ScanForRE(input_str);
                    PrintResult(input_str, res);
                }
            }
            else
            {
                fprintf(stderr, "Input sequence contains invalid entries: %s\n", input_str);
                fprintf(stderr, "Please check the sequence and try again \n");
            }
        }
        else
        {
            fprintf(stderr, "Incorrect Response. Please Enter the Correct Choice\n\n");
        }
    }
}
