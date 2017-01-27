/*
 * visualize_spline_wavelets.cpp (version 1.0)
 * Intended to create plots of synthesis spline wavelets.
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
 * @param x value
 * @param s shift
 */
template <int k>
double Bspline(double x, double s)
{
    if (s <= x && x < (s + k))
    {
        double a = (x - s) / ((s + k - 1) - s);
        double b = ((s + k) - x) / ((s + k) - (s + 1));

        return (a * Bspline<k-1>(x, s) + b * Bspline<k-1>(x, (s + 1)));
    }
    else { return 0.0; }
};

/*
 * @param x value
 * @param s shift
 */
template <>
double Bspline<1>(double x, double s)
{
    if (s <= x && x < (s + 1)) { return 1.0; } else { return 0.0; }
};



/*
 * data structure for spline wavelets
 */
typedef struct {
    int16_t spline_order;      /* B spline order */
    int16_t vanishing_moments; /* vanishing moments */
    const double* mask;
    int mask_start;
    int mask_end;
    int output_start;      /* lower bound of the support */
    int output_end;        /* upper bound of the support */
}psi_t;



/*
 * Visualize spline wavelets
  */
int main(int argc, char** argv)
{
    cout << "Visualize spline wavelets." << endl
         << "Implemented pairs (spline order, vanishing moments): " << endl;

    /*
     * data of implemented spline wavelets
     */
    psi_t psi, psi22, psi23, psi24, psi32, psi33, psi34, psi42, psi43, psi44, psi52, psi53, psi54, psi62, psi63, psi64;

    // <editor-fold defaultstate="collapsed" desc="data of implemented spline wavelets">
    psi22.spline_order = 2;
    psi22.vanishing_moments = 2;
    const double psi22mask[] = {1./4., 1./2., -3./2., 1./2., 1./4.};
    psi22.mask = psi22mask;
    psi22.mask_start = -1;
    psi22.mask_end = 3;
    psi22.output_start = -2;
    psi22.output_end = 2;

    psi23.spline_order = 2;
    psi23.vanishing_moments = 3;
    const double psi23mask[] = {5./32., 5./16., -45./32., 7./8., 11./32., -3./16., -3./32.};
    psi23.mask = psi23mask;
    psi23.mask_start = -1;
    psi23.mask_end = 5;
    psi23.output_start = -2;
    psi23.output_end = 4;

    psi24.spline_order = 2;
    psi24.vanishing_moments = 4;
    const double psi24mask[] = {-3./64., -3./32., 1./4., 19./32., -45./32., 19./32., 1./4., -3./32., -3./64.};
    psi24.mask = psi24mask;
    psi24.mask_start = -1;
    psi24.mask_end = 7;
    psi24.output_start = -2;
    psi24.output_end = 5;

    psi32.spline_order = 3;
    psi32.vanishing_moments = 2;
    const double psi32mask[] = {5./16., 15./16., -15./16., -1./8., 9./8., 3./16.};
    psi32.mask = psi32mask;
    psi32.mask_start = -1;
    psi32.mask_end = 4;
    psi32.output_start = -2;
    psi32.output_end = 4;

    psi33.spline_order = 3;
    psi33.vanishing_moments = 3;
    const double psi33mask[] = {-3./32., -9./32., 7./32., 45./32., -45./32., -7./32., 9./32., 3./32.};
    psi33.mask = psi33mask;
    psi33.mask_start = -1;
    psi33.mask_end = 6;
    psi33.output_start = -2;
    psi33.output_end = 5;

    psi34.spline_order = 3;
    psi34.vanishing_moments = 4;
    const double psi34mask[] = {-7./128., -21./128., 7./32., 35./32., -105./64., 1./64., 19./32., 3./32., -15./128., -5./128.};
    psi34.mask = psi34mask;
    psi34.mask_start = -1;
    psi34.mask_end = 8;
    psi34.output_start = -2;
    psi34.output_end = 6;

    psi42.spline_order = 4;
    psi42.vanishing_moments = 2;
    const double psi42mask[] = {3./16., 3./4., 5./16., -5./2., 5./16., 3./4., 3./16.};
    psi42.mask = psi42mask;
    psi42.mask_start = -2;
    psi42.mask_end = 4;
    psi42.output_start = -3;
    psi42.output_end = 4;

    psi43.spline_order = 4;
    psi43.vanishing_moments = 3;
    const double psi43mask[] = {7./64., 7./16., 0., -35./16., 35./32., 17./16., -1./8., -5./16., -5./64.};
    psi43.mask = psi43mask;
    psi43.mask_start = -2;
    psi43.mask_end = 6;
    psi43.output_start = -3;
    psi43.output_end = 5;

    psi44.spline_order = 4;
    psi44.vanishing_moments = 4;
    const double psi44mask[] = {-5./128., -5./32., -1./128., 3./4., 35./64., -35./16., 35./64., 3./4., -1./128., -5./32., -5./128.};
    psi44.mask = psi44mask;
    psi44.mask_start = -2;
    psi44.mask_end = 8;
    psi44.output_start = -3;
    psi44.output_end = 6;

    psi52.spline_order = 5;
    psi52.vanishing_moments = 2;
    const double psi52mask[] = {7./32., 35./32., 32./32., -105./32., -35./32., 33./32., 25./32., 5./32.};
    psi52.mask = psi52mask;
    psi52.mask_start = -2;
    psi52.mask_end = 5;
    psi52.output_start = -3;
    psi52.output_end = 4;


    psi53.spline_order = 5;
    psi53.vanishing_moments = 3;
    const double psi53mask[] = {-5./64., -25./64., -13./32., 35./32., 35./16., -35./16., -35./32., 13./32., 25./64., 5./64.};
    psi53.mask = psi53mask;
    psi53.mask_start = -2;
    psi53.mask_end = 7;
    psi53.output_start = -3;
    psi53.output_end = 6;


    psi54.spline_order = 5;
    psi54.vanishing_moments = 4;
    const double psi54mask[] = {-45./1024., -225./1024., -171./1024., 945./1024., 735./512., -1365./512., -315./512., 593./512., 575./1024., -165./1024., -175./1024., -35./1024};
    psi54.mask = psi54mask;
    psi54.mask_start = -2;
    psi54.mask_end = 9;
    psi54.output_start = -3;
    psi54.output_end = 6;

    psi62.spline_order = 6;
    psi62.vanishing_moments = 2;
    const double psi62mask[] = {5./32., 15./16., 7./4., -7./16., -77./16., -7./16., 7./4., 15./16., 5./32.};
    psi62.mask = psi62mask;
    psi62.mask_start = -3;
    psi62.mask_end = 5;
    psi62.output_start = -4;
    psi62.output_end = 4;

    psi63.spline_order = 6;
    psi63.vanishing_moments = 3;
    const double psi63mask[] = {45./512., 135./256., 441./512., -63./64., -987./256., 189./128., 693./256., 25./64., -375./512., -105./256., -35./512.};
    psi63.mask = psi63mask;
    psi63.mask_start = -3;
    psi63.mask_end = 7;
    psi63.output_start = -4;
    psi63.output_end = 6;


    psi64.spline_order = 6;
    psi64.vanishing_moments = 4;
    const double psi64mask[] = {-35./1024., -105./512., -165./512., 235./512., 1827./1024., 63./256., -987./256., 63./256., 1827./1024., 235./512., -165./512., -105./512, -35./1024.};
    psi64.mask = psi64mask;
    psi64.mask_start = -3;
    psi64.mask_end = 9;
    psi64.output_start = -4;
    psi64.output_end = 7;
    // </editor-fold>



    /*
     *  user input of spline order and vanishing moments
     */
    int16_t spline_order;
    int16_t vanishing_moments;

    /* check if user passed spline order and vanishing moments as program parameter */
    if (argc == 3)
    {
        spline_order = atoi(argv[1]);
        vanishing_moments = atoi(argv[2]);
        cout << "\nSpline order: " << spline_order << endl
             << "Vanishing moments: " << vanishing_moments << endl;
    }
    else
    {

        cout << "\nInput spline order {2,3,4,5,6}: ";
        cin >> spline_order;
        cout << "Input vanishing moments {2,3,4}: ";
        cin >> vanishing_moments;
    }

    /*
     * select spline wavelet according to user input
     */
    int identifier = (spline_order << 16) | vanishing_moments;
#if 0
    for (int i = 2; i < 5; i++)
    {
        for (int j = 2; j < 5; j ++)
        {
            cout << "(" << i << ", " << j << "): identifier = " << ( (i << 16)|j ) << endl;
        }
    }
#endif

    switch (identifier)
    {
        case 131074:
            psi = psi22;
            break;
        case 131075:
            psi = psi23;
            break;
        case 131076:
            psi = psi24;
            break;
        case 196610:
            psi = psi32;
            break;
        case 196611:
            psi = psi33;
            break;
        case 196612:
            psi = psi34;
            break;
        case 262146:
            psi = psi42;
            break;
        case 262147:
            psi = psi43;
            break;
        case 262148:
            psi = psi44;
            break;
        case 327682:
            psi = psi52;
            break;
        case 327683:
            psi = psi53;
            break;
        case 327684:
            psi = psi54;
            break;
        case 393218:
            psi = psi62;
            break;
        case 393219:
            psi = psi63;
            break;
        case 393220:
            psi = psi64;
            break;
        default:
            cout << "\nCombination of spline order " << spline_order
                 << " and vanishing moments " << vanishing_moments
                 << " not implemented.\nProgram end." << endl;
            return 1;
            break;
    }

    //cout << *(psi.mask + 2) << endl; /* test output third mask entry */

    /* program parameter */
    const int R = 200; /* resolution, e.g., R = 200 */
    double step_size = ( (double) (psi.mask_end - psi.mask_start) / (double) R );
    double scaling_wavelet[R+1] = {0.0}; /* data container */
    double x_values[R+1] = {0.0};
    int stop_write_index = R;
    double value;


    /* setup output file */
    char filename[250];
    sprintf(filename, "spline_wavelet_%d_%d.m", psi.spline_order, psi.vanishing_moments);
    ofstream ofs;

    /* ask user if output should be written to file */
    char answer;
    bool output_file_open = false;
    cout << "\nWrite output to file '" << filename << "' [y,N]? ";
    cin >> answer;
    if (answer == 'y' || answer == 'Y' )
    {
        ofs.open(filename);
        output_file_open = true;
    }
    else
    {
        ofs.close(); /* just to make sure */
    }

    cout << "\nOutput interval: [" << psi.output_start << ", " << psi.output_end << ")" << endl
         << "Resolution: " << step_size << endl
         << "\nOutput:" << endl;

    /* compute and write x values */
    ofs << "X = [";
    for (int i = 0; i <= R; i++)
    {
        x_values[i] = (double)(psi.output_start) + i * step_size;
        //cout << i << "  " << x_values[i] << endl;

        /* write x values to file */
        if (x_values[i] <= psi.output_end)
        {
            ofs << x_values[i] << " ";
            stop_write_index = i;
        }
    }
    ofs << "];\n\nY = [";


    for (int k = psi.mask_start; k <= psi.mask_end; k++)
    {
        for (int i = 0; i <= R; i++)
        {
            //scaling_wavelet[i] += *(psi.mask + (k-psi.mask_start)) * Bspline<spline_order>(2*x_values[i], -2.0 + k);

            switch (psi.spline_order)
            {
                case 2:
                    value = *(psi.mask + (k-psi.mask_start)) * Bspline<2>(2*x_values[i], -2.0 + k);
                    break;
                case 3:
                    value = *(psi.mask + (k-psi.mask_start)) * Bspline<3>(2*x_values[i], -2.0 + k);
                    break;
                case 4:
                    value = *(psi.mask + (k-psi.mask_start)) * Bspline<4>(2*x_values[i], -2.0 + k);
                    break;
                case 5:
                    value = *(psi.mask + (k-psi.mask_start)) * Bspline<5>(2*x_values[i], -2.0 + k);
                    break;
                case 6:
                    value = *(psi.mask + (k-psi.mask_start)) * Bspline<6>(2*x_values[i], -2.0 + k);
                    break;
                default:
                    cout << "\nSpline order " << psi.spline_order << " not implemented. Program end." << endl;
                    return 2;
                    break;
            }
            scaling_wavelet[i] += value;
        }
    }

    /* output results */
    for (int i = 0; i <= stop_write_index; i++)
    {
        cout << "(" << x_values[i] << ", " << scaling_wavelet[i] << ")" << endl;
        ofs << scaling_wavelet[i] << " ";
    }

    ofs << "];\n\nfigure;\nplot(X,Y, 'b', 'LineWidth', 2);" << endl
        << "axis tight;\nset(gca, 'FontSize', 20);" << endl
        << "title('N_" << psi.spline_order << ": spline wavelet \\psi_" << psi.spline_order << "^" << psi.vanishing_moments << "');" << endl;

    ofs.close();

    cout << "\nCardinal B-spline: N_" << psi.spline_order << endl
         << "Vanishing moments: " << psi.vanishing_moments << endl
         << "Output interval: [" << psi.output_start << ", " << psi.output_end << ")" << endl
         << "Resolution: " << step_size << endl;

    if (output_file_open)
    {
        cout << "\nOutput written to: " << endl
             << filename << endl;
    }

    cout << "\nProgram end." << endl;
    cout.flush();

    return 0;
}
