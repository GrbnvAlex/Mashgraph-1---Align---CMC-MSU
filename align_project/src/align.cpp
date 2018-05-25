#include "align.h"
#include <string>
#include <cmath>

using std::string;
using std::cout;
using std::endl;
using std::tie;
using std::get;
using std::make_tuple;
using std::tuple;

const int OFFSET = 15;
#define pi 3.14159265358979323846 

uint right_value(double x)
{
    uint res;
    if (x<0) {res=0;}
    else if (x>255) {res=255;}
    else {res = static_cast<int>(x);}

    return res;
}

int cross_correlation(Image I1, Image I2, int k, int l)
{
    int res = 0, m = I1.n_rows, n = I1.n_cols;

    for (int i = OFFSET; i < m-OFFSET; ++i) {
        for (int j = OFFSET; j < n-OFFSET; ++j) {
            res += get<0>(I1(i,j)) * get<0>(I2(i+k,j+l));
        }
    }

    return res;
}

void max_cross_correlation(Image I1, Image I2, int *i_max, int *j_max)
{
    // I use MSE instead
    *i_max = 0; * j_max = 0;
    int max = cross_correlation(I1,I2,0,0), tmp;

    for (int i=-OFFSET; i<=OFFSET; ++i) {
        for (int j=-OFFSET; j<=OFFSET; ++j) {
            if (not (i==0 && j==0)) {
                tmp = cross_correlation(I1, I2, i, j);
                if (tmp > max) {
                    max= tmp;
                    *i_max = i;
                    *j_max = j;
                }
            }
        }
    }
}

int mse(Image I1, Image I2, int k, int l)
{
    int res = 0, m = I1.n_rows, n = I1.n_cols;

    for (int i = OFFSET; i < m-OFFSET; ++i) {
        for (int j = OFFSET; j < n-OFFSET; ++j) {
            if ( (i+k>=0) && (i+k<m) && (j+l>=0) && (j+l<n) ) {
                int tmp = get<0>(I1(i,j)) - get<0>(I2(i+k,j+l));
                res += tmp*tmp;
            }
        }
    }

    res = res / n / m;
    return res;
}


void min_mse(Image I1, Image I2, int *i_min, int *j_min)
{
    *i_min=0; *j_min=0;
    int mse_min = mse(I1,I2,0,0), tmp;

    for (int i=-OFFSET; i<=OFFSET; ++i) {
        for (int j=-OFFSET; j<=OFFSET; ++j) {
            tmp = mse(I1, I2, i, j);
            if (tmp < mse_min) {
                mse_min = tmp;
                *i_min = i;
                *j_min = j;
            }
        }
    }
}

void pyramid(Image I1, Image I2, int *i, int *j)
{

    int m = I1.n_rows;
    int n = I1.n_cols;
    int tmpi, tmpj;

    // If the image became too small then we start looking for the shift
    if (m/2 < 400 || n/2 < 400) {
        tmpi = *i; tmpj = *j;
        min_mse(I1, I2, &tmpi, &tmpj);
        *i = tmpi; *j = tmpj;
        return;
    }

    // The image is big enough - we make it smaller, getting a level down
    // in the recursion
    I1 = resize(I1, 0.5);
    I2 = resize(I2, 0.5);

    tmpi = *i; tmpj = *j;
    pyramid(I1, I2, &tmpi, &tmpj);
    *i = tmpi; *j = tmpj;

    *i = *i * 2;
    *j = *j * 2;
    return;
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{

    if (isSubpixel) {
        srcImage = resize(srcImage, subScale);
    }

    int m = srcImage.n_rows / 3;
    int n = srcImage.n_cols;

    Image I1 = srcImage.submatrix(0,0,m,n);
    Image I2 = srcImage.submatrix(m,0,m,n);
    Image I3 = srcImage.submatrix(m*2,0,m,n);
    Image res(m,n);

    int i_I3, j_I3, i_I2, j_I2;

    pyramid(I1, I2, &i_I2, &j_I2);
    pyramid(I1, I3, &i_I3, &j_I3);

    int red, green, blue;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            blue = get<0>(I1(i,j));

            if ( i+i_I2>=0 && i+i_I2<m && j+j_I2>=0 && j+j_I2<n ) {
                green = get<0>(I2(i + i_I2,j + j_I2));
            }
            else {
                green = get<1>(I1(i,j));
            }
            if ( i+i_I3>=0 && i+i_I3<m && j+j_I3>=0 && j+j_I3<n ) {
                red = get<0>(I3(i + i_I3,j + j_I3));
            }
            else {
                red = get<0>(I1(i,j));
            }

            res(i,j) = tie(red,green,blue);
        }
    }

    if (isPostprocessing) {
		if (postprocessingType == "--gray-world") {
			res = gray_world(res);
		}
		else if (postprocessingType == "--unsharp") {
			res = unsharp(res);
		}
		else if (postprocessingType == "--autocontrast") {
			res = autocontrast(res, fraction);
		}
	}

    if (isSubpixel) {
        res = resize(res, 1. / subScale).deep_copy();
    }

    return res;

}


Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {{ -1./6.,  -2./3.,  -1./6.},
                             { -2./3.,  4.+1./3.,  -2./3.},
                             { -1./6.,  -2./3.,  -1./6.}};
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {

    Image res(src_image.n_rows, src_image.n_cols);
    double Sr=0, Sb=0, Sg=0, S;
    int r,g,b;

    for (uint i=0; i<src_image.n_rows; ++i) {
        for (uint j=0; j<src_image.n_cols; ++j) {
            tie(r,g,b) = src_image(i,j);
            Sr += r; Sg += g; Sb += b;
        }
    }

    Sr /= src_image.n_rows * src_image.n_cols;
    Sb /= src_image.n_rows * src_image.n_cols;
    Sg /= src_image.n_rows * src_image.n_cols;

    S = (Sr + Sg + Sb) / 3;
    Sr = S/Sr; Sg = S/Sg; Sb = S/Sb;

    for (uint i=0; i<src_image.n_rows; ++i) {
        for (uint j=0; j<src_image.n_cols; ++j) {
            tie(r,g,b) = src_image(i,j);
            res(i,j) = make_tuple(right_value(r*Sr),right_value(g*Sg),right_value(b*Sb));
        }
    }
    return res;
}

Image resize(Image src_image, double scale) {

    int oldw = src_image.n_cols;
    int oldh = src_image.n_rows;
    int neww = static_cast<int>(src_image.n_cols * scale);
    int newh = static_cast<int>(src_image.n_rows * scale);

    Image res(newh,neww);

    int i, j;
	int h, w;
	float t, u, tmp;
	float d1, d2, d3, d4;
	tuple <uint,uint,uint> p1, p2, p3, p4;

	int r, g, b;

	for (i = 0; i < newh; ++i) {
		tmp = static_cast<float>(i) / (newh - 1) * (oldh - 1);
		h = static_cast<int>(floor(tmp));
		if (h < 0) { h = 0; }
        else if (h >= oldh - 1) { h = oldh - 2; }
		u = tmp - h;

		for (j = 0; j < neww; ++j) {

			tmp = static_cast<float>(j) / (neww - 1) * (oldw - 1);
			w = static_cast<int>(floor(tmp));
			if (w < 0) { w = 0; } 
            else if (w >= oldw - 1) { w = oldw - 2; }
			t = tmp - w;

			d1 = (1 - t) * (1 - u);
			d2 = t * (1 - u);
			d3 = (1 - t) * u;
			d4 =  t * u;

			p1 = src_image(h,w);          // (0,0) 
			p2 = src_image(h + 1, w);     // (1,0)
			p3 = src_image(h, w + 1);     // (0,1) 
			p4 = src_image(h + 1, w + 1); // (1,1)

			r =   static_cast<int>(get<0>(p1) * d1) + static_cast<int>(get<0>(p2) * d2) + static_cast<int>(get<0>(p3) * d3) + static_cast<int>(get<0>(p4) * d4);
			g = static_cast<int>(get<1>(p1) * d1) + static_cast<int>(get<1>(p2) * d2) + static_cast<int>(get<1>(p3) * d3) + static_cast<int>(get<1>(p4) * d4);
			b =  static_cast<int>(get<2>(p1) * d1) + static_cast<int>(get<2>(p2) * d2) + static_cast<int>(get<2>(p3) * d3) + static_cast<int>(get<2>(p4) * d4);

            res(i,j) = make_tuple(r, g, b);
		}
	}

    return res;
}


Image mirror(Image src_image, int radius)
{
    int width = src_image.n_cols;
    int height = src_image.n_rows;
    Image mir (height + 2 * radius, width + 2 * radius);
    for (int i=-radius; i<height + radius; ++i) {
        for (int j=-radius; j<width + radius; ++j) {
            if (i >= 0 && i < height && j >= 0 && j < width) { 
                mir(i+radius,j+radius) = src_image(i,j);
            }
            // left side
            else if (i >= 0      && i < height && j<0)        {mir(i+radius,j+radius) = src_image(i,-j);}
            // bottom side
            else if (i >= height && j >= 0     && j < width)  {mir(i+radius,j+radius) = src_image(2*height-i-1,j);}
            // right side
            else if (i < height  && i >= 0     && j >= width) {mir(i+radius,j+radius) = src_image(i,2*width-j-1);}
            // top side
            else if (i < 0       && j >= 0     && j < width)  {mir(i+radius,j+radius) = src_image(-i,j);}
            // top left corner
            else if (i < 0 && j < 0) {mir(i+radius,j+radius) = src_image(-i,-j);}
            // bottom left corner
            else if (i >= height && j < 0) {mir(i+radius,j+radius) = src_image(2*height-i-1,-j);}
            // bottom right corner
            else if (i >= height && j >= width) {mir(i+radius,j+radius) = src_image(2*height-i-1,2*width-j-1);}
            // top right corner
            else if (i < 0       && j > width) {mir(i+radius,j+radius) = src_image(-i,2*width-j-1);}
        }
    }
    return mir;
}

Image custom(Image src_image, Matrix<double> kernel) 
{
    Image res(src_image.n_rows, src_image.n_cols);
    int k_row = kernel.n_rows, k_col = kernel.n_cols; 
    int radius = std::max(kernel.n_rows/2,kernel.n_cols/2);
    int m = src_image.n_rows, n = src_image.n_cols;
    double sum_r=0, sum_g=0, sum_b=0;
    Image mir = mirror(src_image, radius);

    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            sum_r=0; sum_g=0; sum_b=0;
            int r,g,b;
            for (int k=0; k<k_row; ++k) {
                for (int l=0; l<k_col; ++l) {
                    tie(r,g,b) = mir(i+k,j+l);
                    sum_r += r*kernel(k,l);
                    sum_g += g*kernel(k,l);
                    sum_b += b*kernel(k,l);
                }
            }
                r = right_value(sum_r);
                g = right_value(sum_g);
                b = right_value(sum_b);
                res(i,j) = make_tuple(r,g,b);
        }
    }

    return res;
}

Image autocontrast(Image src_image, double fraction) 
{
    int min = 0, max = 0, r, g, b;
    int h[256] = {0};
    double border = src_image.n_rows * src_image.n_cols * fraction;
    for (uint i=0; i<src_image.n_rows; ++i) {
        for (uint j=0; j<src_image.n_cols; ++j) {
            tie(r,g,b) = src_image(i,j);
            h[static_cast<int>(0.2125*r + 0.7154*g + 0.0721*b)]++;
        }
    }

    int s = 0;
    while (s < border) { s += h[min];     min++; }
    s=0;
    while (s < border) { s += h[255-max]; max++; }
    max = 255 - max;

    Image res(src_image.n_rows, src_image.n_cols);
    for (uint i=0; i<src_image.n_rows; ++i) {
        for (uint j=0; j<src_image.n_cols; ++j) {
            tie(r,g,b) = src_image(i,j);
            int br = round(0.2125*r + 0.7154*g + 0.0721*b);
            if (br<min) {res(i,j) = make_tuple(0,0,0);}
            else if (br>max) {res(i,j) = make_tuple(255,255,255);}
            else {
                r = (r - min) * 255 / (max - min);
                g = (g - min) * 255 / (max - min);
                b = (b - min) * 255 / (max - min);
                res(i,j) = make_tuple(right_value(r), right_value(g), right_value(b));
            }
        }
    }
    return res;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    Matrix<double> kernel(2*radius+1,2*radius+1);
    double sum = 0;

    for (int i=0; i<2*radius+1; ++i) {
        for (int j=0; j<2*radius+1; ++j) {
            double tmp = -(i*i+j*j) / (2*sigma*sigma);
            kernel(i,j) = exp(tmp) / (2*pi*sigma*sigma);
            sum += kernel(i,j);
        }
    }

    for (int i=0; i<2*radius+1; ++i) {
        for (int j=0; j<2*radius+1; ++j) { 
            kernel(i,j) /= sum; 
        }
    }

    Image res = custom(src_image, kernel);
    return res;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    Matrix<double>col(2*radius+1, 1);
    Matrix<double>row(1, 2*radius+1);
    double sum = 0;

    for (int i=0; i<2*radius+1; ++i) {
            double tmp = -i*i / (2*sigma*sigma);
            col(i,0) = exp(tmp) / sqrt(2*pi*sigma*sigma);
            row(0,i) = col(i,0);
            sum += col(i,0);
    }

    for (int i=0; i<2*radius+1; ++i) {
            col(i,0) /= sum; 
            row(0,i) /= sum; 
    }
   
    Image res = custom(custom(src_image, row), col);
    return res;
}

int med_of_array(int mas[], int size)
{
    /* We already know that the array has an odd number
     * of elements. That is why here is described only
     * one way of finding the median.
    */
    for (int i=0; i<size-1; ++i) {
        for (int j=0; j<size-i-1; ++j) {
            if (mas[j] > mas[j+1]) {
                int temp = mas[j];
                mas[j] = mas[j+1];
                mas[j + 1] = temp;
            }
        }
    }
    return (mas[size/2]);
}


Image median(Image src_image, int radius) {
    Image mir = mirror(src_image, radius);
    Image res(src_image.n_rows, src_image.n_cols);
    int mir_height = mir.n_rows, mir_width = mir.n_cols;
    int s = (2*radius+1)*(2*radius+1);
    int red[s], green[s], blue[s];
    int cnt=0;

    for (int i=radius; i<mir_height-radius; ++i) {
        for (int j=radius; j<mir_width-radius; ++j) {

            for (int k=-radius; k<=radius; ++k) {
                for (int l=-radius; l<=radius; ++l) {
                    tie(red[cnt],green[cnt],blue[cnt]) = mir(i+k,j+l);
                    cnt++;
                }
            }
            cnt=0;
            res(i-radius,j-radius) = make_tuple(med_of_array(red,s), med_of_array(green,s), med_of_array(blue,s));
        }
    }

    return res;
}

int med(int hist[256], int median_index)
{
    // Returns the median of the input histogram
    int sum = 0, counter = 0;
    do {
        sum += hist[counter];
        counter++;
    }
    while (sum < median_index && counter < 256);

    return counter-1;
}

Image median_linear(Image src_image, int radius) 
{
    /* X - mirror of the input image;
     * r - radius;
     * H - kernel histogram;
     * Y - an output image;
     * median(H) - median of kernel histogram H;
    */

    Image mir = mirror(src_image, radius);
    Image res(src_image.n_rows, src_image.n_cols);
    int hist_r[256] = {0}, hist_g[256] = {0}, hist_b[256] = {0};
    int med_i = (2*radius+1)*(2*radius+1)/2;
    int r,g,b;
    int mir_height = mir.n_rows, mir_width = mir.n_cols;
    
    for (int i=radius; i<mir_height-radius; ++i) {
        for (int j=radius; j<mir_width-radius; ++j) {

            if (j == radius) {
                // Zero fill histogram H
                for (int k=0; k<256; ++k) {
                    hist_r[k] = 0; hist_g[k] = 0; hist_b[k] = 0;
                }

                // Initialize kernel histogram H
                for (int k=-radius; k<=radius; ++k) {
                    for (int l=-radius; l<=radius; ++l) {
                        tie(r,g,b) = mir(i+k,j+l);
                        hist_r[r]++; hist_g[g]++; hist_b[b]++;
                    }
                }
                // Y[i,j] <- median(H)
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }
            else {
                for (int k=-radius; k<=radius; ++k) {
                    // Remove X[i+k,j-r-1] from H
                    tie(r,g,b) = mir(i+k,j-radius-1);
                    hist_r[r]--; hist_g[g]--; hist_b[b]--;
                    // Add X[i+k,j+r] to H
                    tie(r,g,b) = mir(i+k,j+radius);
                    hist_r[r]++; hist_g[g]++; hist_b[b]++;
                    
                }
                // Y[i,j] <- median(H)
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }
        }
    }
    return res;
}

Image median_const(Image src_image, int radius) {
    Image mir = mirror(src_image, radius);
    Image res(src_image.n_rows, src_image.n_cols);
    int rows = mir.n_rows, cols = mir.n_cols;
    int hist_r[256] = {0}, hist_g[256] = {0}, hist_b[256] = {0};
    int h_r[mir.n_cols][256];
    int h_g[mir.n_cols][256];
    int h_b[mir.n_cols][256];

    for (uint i=0; i<mir.n_cols; i++) {
        for (uint j=0; j<255; j++) {
            h_r[i][j] = 0;
            h_g[i][j] = 0;
            h_b[i][j] = 0;
        }
    }

    int med_i = (2*radius+1)*(2*radius+1)/2;
    int r,g,b;

    for (int i=radius; i<rows-radius; ++i) {
        for (int j=radius; j<cols-radius; ++j) {

            // 1 row; 1 col
            if (i == radius && j == radius) {

                // 1) Filling (2*radius+1) first column histograms;
                for (int k=-radius; k<=radius; ++k) {
                    for (int l=-radius; l<=radius; ++l) {
                        tie(r,g,b) = mir(i+k,j+l);
                        h_r[j+l][r]++; h_g[j+l][g]++; h_b[j+l][b]++;
                    }
                }

                // 2) hist = hist + (2*radius+1) first column histograms;
                for (int k=0; k<2*radius+1; ++k) {
                    for (int l=0; l<256; ++l) {
                        hist_r[l] += h_r[k][l]; hist_g[l] += h_g[k][l]; hist_b[l] += h_b[k][l];
                    }
                }

                // 3) Forming the res(i,j) using hist;
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }

            // 1 row; 2 col ... (src_image.n_cols-1) col
            else if (i == radius) {

                // 1) Filling the next column histogram;
                for (int k=-radius; k<=radius; ++k) {
                    tie(r,g,b) = mir(i+k,j+radius);
                    h_r[j+radius][r]++; h_g[j+radius][g]++; h_b[j+radius][b]++;
                }

                // 2) hist = hist + that filled column histogram - previous;
                for (int k=-radius; k<=radius; ++k) {
                    for (int l=0; l<256; ++l) {
                        hist_r[l] = hist_r[l] + h_r[j+radius][l] - h_r[j-radius-1][l];
                        hist_g[l] = hist_g[l] + h_g[j+radius][l] - h_g[j-radius-1][l];
                        hist_b[l] = hist_b[l] + h_b[j+radius][l] - h_b[j-radius-1][l];
                    } 
                }

                // 3) Forming the res(i,j) using hist;
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }

            // 2 row ... (src_image.n_rows-1); 1 col
            else if (j == radius) {

                // 1) Moving on a row down (2*radius+1) first column histograms;
                for (int l=-radius; l<=radius; ++l) {
                    // Minus an old pixel;
                    tie(r,g,b) = mir(i-radius-1,j+l);
                    h_r[j+l][r]--; h_g[j+l][g]--; h_b[j+l][b]--;
                    // Plus a new pixel;
                    tie(r,g,b) = mir(i+radius,j+l);
                    h_r[j+l][r]++; h_g[j+l][g]++; h_b[j+l][b]++;
                }

                // 2) Zero fill hist
                for (int k=0; k<256; ++k) {
                    hist_r[k] = 0; hist_g[k] = 0; hist_b[k] = 0;
                }

                // 3) hist = hist + (2*radius+1) first column histograms;
                for (int k=0; k<2*radius+1; ++k) {
                    for (int l=0; l<256; ++l) {
                        hist_r[l] += h_r[k][l]; hist_g[l] += h_g[k][l]; hist_b[l] += h_b[k][l];
                    }
                }

                // 4) Forming the res(i,j) using hist;
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }

            // 2 row ... (src_image.n_rows-1); 2 col ... (src_image.n_cols-1) col
            else {

                // 1) Moving on a row down the next column histogram;
                // Minus an old pixel;
                tie(r,g,b) = mir(i-radius-1,j+radius);
                h_r[j+radius][r]--; h_g[j+radius][g]--; h_b[j+radius][b]--;
                // Plus a new pixel;
                tie(r,g,b) = mir(i+radius,j+radius);
                h_r[j+radius][r]++; h_g[j+radius][g]++; h_b[j+radius][b]++;

                // 2) hist = hist + that new moved down column histogram - previous;
                for (int k=-radius; k<=radius; ++k) {
                    for (int l=0; l<256; ++l) {
                        hist_r[l] = hist_r[l] + h_r[j+radius][l] - h_r[j-radius-1][l];
                        hist_g[l] = hist_g[l] + h_g[j+radius][l] - h_g[j-radius-1][l];
                        hist_b[l] = hist_b[l] + h_b[j+radius][l] - h_b[j-radius-1][l];
                    } 
                }

                // 3) Forming the res(i,j) using hist;
                res(i-radius,j-radius) = make_tuple(med(hist_r,med_i), med(hist_g,med_i), med(hist_b,med_i));
            }
        
        }
    }

    return res;
}

Image abs_grad(Image one, Image two)
{
    int r,g,b;
    Image res(one.n_rows, one.n_cols);
    for (uint i=0; i<one.n_rows; ++i) {
        for (uint j=0; j<one.n_cols; ++j) {
            r = sqrt(  get<0>(one(i,j))*get<0>(one(i,j)) + get<0>(two(i,j))*get<0>(two(i,j))  );
            g = sqrt(  get<1>(one(i,j))*get<1>(one(i,j)) + get<1>(two(i,j))*get<1>(two(i,j))  );
            b = sqrt(  get<2>(one(i,j))*get<0>(one(i,j)) + get<2>(two(i,j))*get<2>(two(i,j))  );
            res(i,j) = make_tuple(right_value(r),right_value(g),right_value(b));
        }
    }
    return res;
}

void choose_neighbour(double x, int *i, int *j)
{
    // x in [-pi,pi] (from atan2)

    if      (x >= -pi/8 && x < pi/8.)           {*i=0; *j=1;}    // right
    else if (x >= pi/8. && x < 3*pi/8.)       {*i=-1; *j=1;}   // top right corner
    else if (x >= 3*pi/8. && x < 5*pi/8.)     {*i=-1; *j=0;}   // top
    else if (x >= 5*pi/8. && x < 7*pi/8.)        {*i=-1; *j=-1;}  // top left corner
    else if (x >= 3*pi/8. && x < -pi/8.)          {*i=1; *j=1;}    // bot right corner
    else if (x >= -5*pi/8. && x < -3*pi/8.)     {*i=1; *j=0;}    // bot
    else if (x >= -7*pi/8. && x < -5*pi/8.)   {*i=1; *j=-1;}   // bot left corner
    else                                    {*i=0; *j=-1;}   // left


}

Image supression_nonmax(Image grad, Image one, Image two)
{
    Image mir = mirror(grad,1);
    Image res(grad.n_rows, grad.n_cols);
    double grad_r,grad_g,grad_b;
    int r,g,b;
    int k,l;
    for (uint i=0; i<grad.n_rows; ++i) {
        for (uint j=0; j<grad.n_cols; ++j) {
            grad_r = atan2(  get<0>(two(i,j)), get<0>(one(i,j))  );
            grad_g = atan2(  get<1>(two(i,j)), get<1>(one(i,j))  );
            grad_b = atan2(  get<2>(two(i,j)), get<2>(one(i,j))  );

            choose_neighbour(grad_r, &k, &l); 
            // +1 because we deal with mirror
            if (get<0>(mir(i+1,j+1)) > get<0>(mir(i+k+1,j+l+1)) && get<0>(mir(i+1,j+1)) > get<0>(mir(i-k+1,j-l+1)) ) { 
                r = get<0>(mir(i+1,j+1));
            } else { r = 0; g=0; b=0; }

            choose_neighbour(grad_g, &k, &l); 
            if (get<1>(mir(i+1,j+1)) > get<1>(mir(i+k+1,j+l+1)) && get<1>(mir(i+1,j+1)) > get<1>(mir(i-k+1,j-l+1)) ) { 
                g = get<1>(mir(i+1,j+1)); 
            } else { g = 0; }

            choose_neighbour(grad_b, &k, &l); 
            if (get<2>(mir(i+1,j+1)) > get<2>(mir(i+k+1,j+l+1)) && get<2>(mir(i+1,j+1)) > get<2>(mir(i-k+1,j-l+1)) ) { 
                b = get<2>(mir(i+1,j+1)); 
            } else { b = 0; }

            res(i,j) = make_tuple(r,g,b);
        }
    }
    return res;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    // 1 - Smoothing the image
    Image res = gaussian(src_image, 1.4, 2);

    // 2 - Computing gradients
    Image Ix = sobel_x(res);
    Image Iy = sobel_y(res);
    Image Gr = abs_grad(Ix, Iy);
    Image Sn = supression_nonmax(Gr, Ix, Iy);

    // He finds the absolute mean of the image normally
    // But he doesn't sharpen the border - most of the border disappears after nonmax supression =(
    // Although the algotithm looks pretty smooth

    return Sn;
}