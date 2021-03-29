using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Imaging;
using Microsoft.Win32;
using System.Drawing;
using System.IO;
using System.Drawing.Imaging;
using System.Windows.Input;
using System.Text.RegularExpressions;

namespace ImageFilterApp
{
    public partial class MainWindow : Window
    {

        private BitmapImage original_Image;
        private BitmapImage previous_Image;

        public MainWindow()
        {
            InitializeComponent();
        }

        // The following 5 functions were taken from docs.microsoft or stackoverflow to help with some basic functionalities

        private void Load_Image(object sender, RoutedEventArgs e)
        {
            OpenFileDialog op = new OpenFileDialog();
            op.Title = "Select a picture";
            op.Filter = "All supported graphics|    *.jpg;*.jpeg;*.png|" +
              "JPEG (*.jpg;*.jpeg)|*.jpg;*.jpeg|" +
              "Portable Network Graphic (*.png)|*.png";
            if (op.ShowDialog() == true)
            {
                input_Image.Source = new BitmapImage(new Uri(op.FileName));
                original_Image = new BitmapImage(new Uri(op.FileName));
            }

        }

        private void Save_Image(object sender, RoutedEventArgs e)
        {
            SaveFileDialog dlg = new SaveFileDialog();
            dlg.FileName = "Untitled";
            dlg.Filter = "JPeg Image|*.jpg|Bitmap Image|*.bmp";
            if (dlg.ShowDialog() == true)
            {
                var encoder = new JpegBitmapEncoder(); // Or PngBitmapEncoder, or whichever encoder you want
                encoder.Frames.Add(BitmapFrame.Create(input_Image.Source as BitmapImage));
                using (var stream = dlg.OpenFile())
                {
                    encoder.Save(stream);
                }
            }
        }

        private Bitmap BitmapImage2Bitmap(BitmapImage bitmapImage)
        {
            // BitmapImage bitmapImage = new BitmapImage(new Uri("../Images/test.png", UriKind.Relative));

            using (MemoryStream outStream = new MemoryStream())
            {
                BitmapEncoder enc = new BmpBitmapEncoder();
                enc.Frames.Add(BitmapFrame.Create(bitmapImage));
                enc.Save(outStream);
                System.Drawing.Bitmap bitmap = new System.Drawing.Bitmap(outStream);

                return new Bitmap(bitmap);
            }
        }

        private BitmapImage ToBitmapImage(Bitmap bitmap)
        {
            using (var memory = new MemoryStream())
            {
                bitmap.Save(memory, ImageFormat.Png);
                memory.Position = 0;

                var bitmapImage = new BitmapImage();
                bitmapImage.BeginInit();
                bitmapImage.StreamSource = memory;
                bitmapImage.CacheOption = BitmapCacheOption.OnLoad;
                bitmapImage.EndInit();
                bitmapImage.Freeze();

                return bitmapImage;
            }
        }

        private static T Clamp<T>(T val, T min, T max) where T : IComparable<T>
        {
            if (val.CompareTo(min) < 0) return min;
            else if (val.CompareTo(max) > 0) return max;
            else return val;
        }

        //===============================================================================


        // Function filters
        private void Inversionfilter()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);
                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, 255 - pixelColor.R, 255 - pixelColor.G, 255 - pixelColor.B));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void BrightnesCorrection()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);

            int correction = 30;
            int red, green, blue;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);

                    red = Clamp(pixelColor.R + correction, 0, 255);
                    green = Clamp(pixelColor.G + correction, 0, 255);
                    blue = Clamp(pixelColor.B + correction, 0, 255);

                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, red, green, blue));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void ContrastEnhancement()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            double aCoeff = 1.65, bCoeff = -82.46;
            int red, green, blue;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);

                    red = Clamp(Convert.ToInt32(aCoeff * pixelColor.R + bCoeff), 0, 255);
                    green = Clamp(Convert.ToInt32(aCoeff * pixelColor.G + bCoeff), 0, 255);
                    blue = Clamp(Convert.ToInt32(aCoeff * pixelColor.B + bCoeff), 0, 255);

                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, red, green, blue));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private void GammaCorrection()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            double gamma = 1.5;

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    Color pixelColor = pic.GetPixel(x, y);
                    pic.SetPixel(x, y, Color.FromArgb(pixelColor.A, 
                        GammaCorrectionCalculation(gamma,pixelColor.R), 
                        GammaCorrectionCalculation(gamma, pixelColor.G),
                        GammaCorrectionCalculation(gamma, pixelColor.B)));
                }
            }
            output_Image.Source = ToBitmapImage(pic);
        }

        private int GammaCorrectionCalculation(double gamma, double color)
        {
            return Convert.ToInt32(Math.Pow(color / 255, gamma) * 255);
        }

        // ===============================================================================

        // Convolution filters

        private void Blur3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelCount = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for(int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width) 
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += pixel.R;
                                avgG += pixel.G;
                                avgB += pixel.B;

                                blurPixelCount++;
                            }
                        }
                    }

                    avgR = avgR / blurPixelCount;
                    avgG = avgG / blurPixelCount;
                    avgB = avgB / blurPixelCount;

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void GaussianSmoothing3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { 0, 1, 0 }, 
                                                  { 1, 4, 1 },
                                                  { 0, 1, 0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = avgR / blurPixelWeighs;
                    avgG = avgG / blurPixelWeighs;
                    avgB = avgB / blurPixelWeighs;

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void Sharpen3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);

            int b = 5, a = 1;
            int s = b - 4*a;

            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            double[,] kernelWeighs = new double[3, 3] { {   0, -a/s,   0 },
                                                        { -a/s, b/s, -a/s },
                                                        {   0, -a/s,   0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    double blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                Clamp(avgR += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R), 0, 255);
                                Clamp(avgG += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G), 0, 255);
                                Clamp(avgB += Convert.ToInt32(kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B), 0, 255);

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                        }
                            }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / blurPixelWeighs),0,255);
                    avgG = Clamp(Convert.ToInt32(avgG / blurPixelWeighs),0,255);
                    avgB = Clamp(Convert.ToInt32(avgB / blurPixelWeighs),0,255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void EdgeDetect3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { 0, -1,  0 },
                                                  { 0,  1,  0 },
                                                  { 0,  0,  0 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / 1), 0, 255);
                    avgG = Clamp(Convert.ToInt32(avgG / 1), 0, 255);
                    avgB = Clamp(Convert.ToInt32(avgB / 1), 0, 255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        private void Emboss3x3()
        {
            Bitmap pic = BitmapImage2Bitmap(input_Image.Source as BitmapImage);
            Bitmap output = new Bitmap(pic.Width, pic.Height);
            int kernelSize = 3;
            int matrixLimits = Convert.ToInt32((kernelSize - 1) / 2);
            int[,] kernelWeighs = new int[3, 3] { { -1, 0, 1 },
                                                  { -1, 1, 1 },
                                                  { -1, 0, 1 } };

            for (int y = 0; (y <= (pic.Height - 1)); y++)
            {
                for (int x = 0; (x <= (pic.Width - 1)); x++)
                {
                    int avgR = 0, avgG = 0, avgB = 0;
                    int blurPixelWeighs = 0;

                    int kernelCoeffX = 0, kernelCoeffY = 0;

                    for (int ky = -matrixLimits; ky <= matrixLimits; ky++)
                    {
                        for (int kx = -matrixLimits; kx <= matrixLimits; kx++)
                        {
                            if (y + ky >= 0 && y + ky < pic.Height && x + kx >= 0 && x + kx < pic.Width)
                            {
                                Color pixel = pic.GetPixel(x + kx, y + ky);

                                avgR += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.R;
                                avgG += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.G;
                                avgB += kernelWeighs[kernelCoeffX, kernelCoeffY] * pixel.B;

                                blurPixelWeighs += kernelWeighs[kernelCoeffX, kernelCoeffY];
                                kernelCoeffX++;
                            }
                        }
                        kernelCoeffX = 0;
                        kernelCoeffY++;
                    }

                    avgR = Clamp(Convert.ToInt32(avgR / 1), 0, 255);
                    avgG = Clamp(Convert.ToInt32(avgG / 1), 0, 255);
                    avgB = Clamp(Convert.ToInt32(avgB / 1), 0, 255);

                    output.SetPixel(x, y, Color.FromArgb(pic.GetPixel(x, y).A, avgR, avgG, avgB));
                }
            }

            output_Image.Source = ToBitmapImage(output);
        }

        // ===============================================================================

        // Button events 

        private void Apply_Function(object sender, RoutedEventArgs e)
        {
            previous_Image = input_Image.Source as BitmapImage;
            input_Image.Source = output_Image.Source;
        }

        private void Undo_Click(object sender, RoutedEventArgs e)
        {
            input_Image.Source = previous_Image;
        }

        private void UndoAll_Click(object sender, RoutedEventArgs e)
        {
            input_Image.Source = original_Image;
        }

        private void Inversion_Click(object sender, RoutedEventArgs e)
        {
            Inversionfilter();
        }

        private void Brightnes_Click(object sender, RoutedEventArgs e)
        {
            BrightnesCorrection();
        }

        private void Contrast_Click(object sender, RoutedEventArgs e)
        {
            ContrastEnhancement();
        }

        private void Gamma_Click(object sender, RoutedEventArgs e)
        {
            GammaCorrection();
        }

        private void Blur_Click(object sender, RoutedEventArgs e)
        {
            Blur3x3();
        }

        private void GaussianBlur_Click(object sender, RoutedEventArgs e)
        {
            GaussianSmoothing3x3();
        }

        private void Sharpen_Click(object sender, RoutedEventArgs e)
        {
            Sharpen3x3();
        }

        private void Edge_Click(object sender, RoutedEventArgs e)
        {
            EdgeDetect3x3();
        }

        private void Emboss_Click(object sender, RoutedEventArgs e)
        {
            Emboss3x3();
        }

        // =============================================================

    }
}
