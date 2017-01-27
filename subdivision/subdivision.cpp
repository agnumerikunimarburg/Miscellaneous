/*
 * subdivision.cpp (version 1.0)
 * Intended to create plots of refinable functions.
 *
 * This software is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either expressed or implied.
 *
 * Contact:  AG Numerik, Philipps-University Marburg
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdint.h>

using namespace std;



/*
 * data structure for refinable mask
 */
typedef struct {
    const double* entry;
    int length;
    string name;
}mask_t;



/*
 *
 */
int main(int argc, char** argv)
{
    cout << "Subdivision to visualize refinable functions." << endl;

    /*
     * data of implemented refinable functions
     */
    const int mask_max_length = 20; /* maximal length of mask */
    mask_t mask, maskHaar, maskN2, maskN3, maskN4, maskN5, maskN6, maskI4, maskI6, maskD2, maskCDF13d, maskCDF22d, maskCDF35d, maskCDF46d;

    // <editor-fold defaultstate="collapsed" desc="data of implemented refinable functions">
    const double Haar[] = {1.0, 1.0}; /* Haar */
    maskHaar.entry = Haar;
    maskHaar.length = 2;
    maskHaar.name = "Haar";

    const double N2[] = {1.0/2.0, 1.0, 1.0/2.0}; /* N_2 */
    maskN2.entry = N2;
    maskN2.length = 3;
    maskN2.name = "N_2";

    const double N3[] = {1.0/4.0, 3.0/4.0, 3.0/4.0, 1.0/4.0}; /* N_3 */
    maskN3.entry = N3;
    maskN3.length = 4;
    maskN3.name = "N_3";

    const double N4[] = {1.0/8.0, 1.0/2.0, 3.0/4.0, 1.0/2.0, 1.0/8.0}; /* N_4 */
    maskN4.entry = N4;
    maskN4.length = 5;
    maskN4.name = "N_4";

    const double N5[] = {1.0/16.0, 5.0/16.0, 10.0/16.0, 10.0/16.0, 5.0/16.0, 1.0/16.0}; /* N_5 */
    maskN5.entry = N5;
    maskN5.length = 6;
    maskN5.name = "N_5";

    const double N6[] = {1.0/32.0, 6.0/32.0, 15.0/32.0, 20.0/32.0, 15.0/32.0, 6.0/32.0, 1.0/32.0}; /* N_6 */
    maskN6.entry = N6;
    maskN6.length = 7;
    maskN6.name = "N_6";

    const double I4[] = {-1.0/16.0, 0.0, 9.0/16.0, 1.0, 9.0/16.0, 0.0, -1.0/16.0}; /* I_4 */
    maskI4.entry = I4;
    maskI4.length = 7;
    maskI4.name = "I_4";

    const double I6[] = {3.0/256.0, 0.0, -25.0/256.0, 0.0, 150.0/256.0, 1.0, 150.0/256.0, 0.0, -25.0/256.0, 0.0, 3.0/256.0}; /* I_6 */
    maskI6.entry = I6;
    maskI6.length = 11;
    maskI6.name = "I_6";

    const double D2[] = {(1.0+sqrt(3.0))/4.0, (3.0+sqrt(3.0))/4.0, (3.0-sqrt(3.0))/4.0, (1.0-sqrt(3.0))/4.0}; /* Daubechies 2 */
    maskD2.entry = D2;
    maskD2.length = 4;
    maskD2.name = "Daubechies_2";

    const double CDF13d[] = {-1.0/8.0, 1.0/8.0, 1.0, 1.0, 1.0/8.0, -1.0/8.0}; /* CDF (1,3) dual */
    maskCDF13d.entry = CDF13d;
    maskCDF13d.length = 6;
    maskCDF13d.name = "CDF_1_3_dual";

    const double CDF22d[] = {-1.0/4.0, 1.0/2.0, 3.0/2.0, 1.0/2.0, -1.0/4.0}; /* CDF (2,2) dual */
    maskCDF22d.entry = CDF22d;
    maskCDF22d.length = 5;
    maskCDF22d.name = "CDF_2_2_dual";

    const double CDF35d[] = {-5.0/256.0, 15.0/256.0, 19.0/256.0, -97.0/256.0, -26.0/256.0, 350.0/256.0, 350.0/256.0, -26.0/256.0, -97.0/256.0, 19.0/256.0, 15.0/256.0, -5.0/256.0}; /* CDF (3,5) dual */
    maskCDF35d.entry = CDF35d;
    maskCDF35d.length = 12;
    maskCDF35d.name = "CDF_3_5_dual";

    const double CDF46d[] = {70.0/8192.0, -70.0/2048.0, -110.0/8192.0, 230.0/1024.0, -1114.0/8192.0, -1466.0/2048.0, 5250.0/8192.0, 1050.0/512.0, 5250.0/8192.0, -1466.0/2048.0, -1114.0/8192.0, 230.0/1024.0, -110.0/8192.0, -70.0/2048.0, 70.0/8192.0}; /* CDF (4,6) dual */
    maskCDF46d.entry = CDF46d;
    maskCDF46d.length = 15;
    maskCDF46d.name = "CDF_4_6_dual";
    // </editor-fold>

    int selected_mask = 0;

    cout << "\nNumber | Refinement function" << endl
         << "----------------------------" << endl
         << "   1   |  Haar            " << endl
         << "   2   |  N_2             " << endl
         << "   3   |  N_3             " << endl
         << "   4   |  N_4             " << endl
         << "   5   |  N_5             " << endl
         << "   6   |  N_6             " << endl
         << "   7   |  I_4             " << endl
         << "   8   |  I_6             " << endl
         << "   9   |  Daubechies 2    " << endl
         << "  10   |  CDF (1,3) dual  " << endl
         << "  11   |  CDF (2,2) dual  " << endl
         << "  12   |  CDF (3,5) dual  " << endl
         << "  13   |  CDF (4,6) dual  " << endl
         << "\nSelect mask {1, 2, ...}: ";
    cin >> selected_mask;

    switch (selected_mask)
    {
        // <editor-fold defaultstate="collapsed" desc="switch cases">
        case 1:
            mask = maskHaar;
            break;
        case 2:
            mask = maskN2;
            break;
        case 3:
            mask = maskN3;
            break;
        case 4:
            mask = maskN4;
            break;
        case 5:
            mask = maskN5;
            break;
        case 6:
            mask = maskN6;
            break;
        case 7:
            mask = maskI4;
            break;
        case 8:
            mask = maskI6;
            break;
        case 9:
            mask = maskD2;
            break;
        case 10:
            mask = maskCDF13d;
            break;
        case 11:
            mask = maskCDF22d;
            break;
        case 12:
            mask = maskCDF35d;
            break;
        case 13:
            mask = maskCDF46d;
            break;
        // </editor-fold>

        default:
            cout << "\nNo mask implemented for this input.\nProgram end." << endl;
            return 1;
            break;
    }

    cout << "\nOutput:" << endl;

    /* program parameters */
    const int max_steps = 8; /* maximal subdivision steps, e.g., 8 */
    const int M = (mask.length) - 1; /* length of mask -1 */
    //double S[max_steps + 2][(int) (pow(2.0, (max_steps + 1)) * (double) M) + 1] = {0.0};
    double S[max_steps + 2][(int) (pow(2.0, (max_steps + 1)) * (double) mask_max_length)] = {0.0};
    int lowerk, upperk, lowerm, upperm;
    double value;


    /* subdivision scheme */
    S[1][M+1] = 1.0;

    for (int j = 1; j <= max_steps; j++) /* subdivision steps */
    {
        lowerk = (-1) * (int) (pow(2.0, (j)) * (double) M);
        upperk =        (int) (pow(2.0, (j)) * (double) M);

        for (int k = lowerk; k <= upperk; k++)
        {
            value = 0.0;
            lowerm = (int) max( (-1.0) * pow(2.0, (j-1)) * (double) M, ceil((double) (k-M) / 2.0) );
            upperm = (int) min( pow(2.0, (j-1)) * (double) M, floor((double) k / 2.0) );

            for (int m = lowerm; m <= upperm; m++)
            {
                //value += mask[k - 2*m] * S[j][m + (int) pow(2.0, (j-1)) * M + 1];
                value += *(mask.entry + (k - 2*m)) * S[j][m + (int) pow(2.0, (j-1)) * M + 1];
            }

            S[j+1][k + (int) (pow(2.0, j) * (double) M) + 1] =  value;

            cout << "(" << j << ", " << k << "): value = " << value << endl;
        }
    }


    /* setup output file */
    char filename[250];
    const char *cstr = mask.name.c_str();
    sprintf(filename, "subdivision_%s.m", cstr); //, max_steps);
    ofstream ofs;

    /* ask user if output should be written to file */
    char answer;
    int counter = 0;
    cout << "\nWrite output to file '" << filename << "' [y,N]? ";
    cin >> answer;
    if (answer == 'y' || answer == 'Y' )
    {
        ofs.open(filename);

       // <editor-fold defaultstate="collapsed" desc="results are written to file">

        /* plot of every step */
#if 1
        ofs << "figure;" << endl;

        for (int step = 0; step < max_steps+1; step ++)
        {
            counter = 0;
            ofs << "\nX = [";
            for (double d = 0; d < M; d += pow(2.0, (-1)*step))
            {
                ofs << d << " ";
                counter++;
            }
            ofs << "];\nY = [";
            //for (int j = 1+pow(2.0,max_steps)*M; j < sizeof(S[max_steps+1])/sizeof(*S[max_steps+1]); j++ )
            for (int j = 1+pow(2.0,step)*M; j < (1+pow(2.0,step)*M) + counter; j++ )
            {
                ofs << S[step+1][j] << " ";
            }
            ofs << "];\nplot(X,Y);\ntitle('step " << step << "');\npause(0.4);clf;" << endl;
            //cout << counter << endl;
        }

        ofs << "\nplot(X,Y, 'b', 'LineWidth', 2);" << endl
            << "%axis tight;\n%set(gca, 'FontSize', 14);" << endl
            << "title('mask = ";
        for (int i = 0; i < M; i++)
        {
            //ofs << mask[i] << ", ";
            ofs << *(mask.entry + i) << ", ";
        }
        ofs << *(mask.entry + M) //<< mask[M]
            << "');" << endl
            << "\nset(gca, 'FontSize', 20);" << endl;


        /* plot of result only */
#else
        ofs << "X = [";
        for (double d = 0; d < M; d += pow(2.0, (-1)*max_steps))
        {
            ofs << d << " ";
            counter++;
        }
        ofs << "];\nY = [";
        //for (int j = 1+pow(2.0,max_steps)*M; j < sizeof(S[max_steps+1])/sizeof(*S[max_steps+1]); j++ )
        for (int j = 1+pow(2.0,max_steps)*M; j < (1+pow(2.0,max_steps)*M) + counter; j++ )
        {
            ofs << S[max_steps+1][j] << " ";
        }


        ofs << "];\n\nfigure;\nplot(X,Y); % , 'k', 'LineWidth', 2);" << endl
            << "%axis tight;\n%set(gca, 'FontSize', 14);" << endl
            << "title('mask = ";
        for (int i = 0; i < M; i++)
        {
            //ofs << mask[i] << ", ";
            ofs << *(mask.entry + i) << ", ";
        }
        ofs << *(mask.entry + M) //<< mask[M]
            << "');" << endl;
#endif

        // </editor-fold>

        ofs.close();

        cout << "\nOutput written to: " << endl
             << filename << endl;
    }

    /* check norm of mask: the sum of the mask elements should be equal to 2 */
    double norm_of_mask = 0.0;
    for (int i = 0; i < mask.length; i++)
    {
        norm_of_mask += *(mask.entry + i);
    }

    if (norm_of_mask != 2.0)
    {
        cout << "\nWarning: Norm of mask does not equal 2." << endl;
    }

    cout << "\nProgram end." << endl;
    cout.flush();

    return 0;
}
