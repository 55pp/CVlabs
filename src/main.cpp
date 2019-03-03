#include <QCoreApplication>
#include <iostream>
#include <QBitmap>
#include <QVector>
#include <functional>
#include <algorithm>
#include <cvfunctions.h>

using namespace std::placeholders;

using Float = long double;
using Vektor = QVector<Float>;
using PointsVector = QVector<QPoint>;
using Matrix = QVector<Vektor>;

using PSideFunction = int (*)(const Matrix&, int, int);

struct ImageProcessing
{
  Matrix img;
  Matrix core;
  SideMode sideMode;
  CoreMode coreMode;
};

struct ElementOfOctave
{
  ImageProcessing imgProc;
  int numberOfOctave;
  double localSigma;
  double globalSigma;
};

template <typename T>
int
toByte(T num)
{
  return (num > 255 ? 255 : (num < 0 ? 0 : num));
}

inline void
outstr(const QString& msg)
{
  std::cout << msg.toStdString() << std::endl;
}

inline void
debugOutput(const QString& msg)
{
#ifdef QT_DEBUG
  std::cout << msg.toStdString() << std::endl;
#endif
}

QImage
loadImage(const QString& fileName)
{
  QImage image;
  if (!image.load(fileName)) {
    outstr("File " + fileName + " not open");
  }

  return image;
}

double
rgbToGrayscale(const QRgb& color)
{
  QColor rgb(color);
  return (0.213 * rgb.red() + 0.715 * rgb.green() + 0.072 * rgb.blue());
}

int
bounded(int num, int bound)
{
  int res;

  if (num < 0) {
    res = 0;
  } else if (num >= bound) {
    res = bound - 1;
  } else {
    res = num;
  }

  return res;
}

int
reflected(int num, int bound)
{
  int res;

  if (num < 0) {
    res = -(num + 1);
  } else if (num >= bound) {
    res = bound + bound - 1 - num;
  } else {
    res = num;
  }

  return res;
}

int
rounded(int num, int bound)
{
  int res;

  if (num < 0) {
    res = bound + num;
  } else if (num >= bound) {
    res = num % bound;
  } else {
    res = num;
  }

  return res;
}

int
blackSide(const Matrix& image, int x, int y)
{
  (void)image;
  (void)x;
  (void)y;
  return 0;
}

int
roundSide(const Matrix& image, int x, int y)
{
  int color = image[rounded(x, image.size())][rounded(y, image[0].size())];
  return color;
}

int
boundSide(const Matrix& image, int x, int y)
{
  int color = image[bounded(x, image.size())][bounded(y, image[0].size())];
  return color;
}

int
reflectSide(const Matrix& image, int x, int y)
{
  int color = image[reflected(x, image.size())][reflected(y, image[0].size())];
  return color;
}

Matrix
rgbImageToGrayscaleImage(const QImage& image)
{
  Matrix grayImage;
  grayImage.fill(*new Vektor(image.width()), image.height());
  for (int i(0); i < image.height(); ++i) {
    for (int j(0); j < image.width(); ++j) {
      QRgb color = image.pixel(j, i);
      double sum = rgbToGrayscale(color);
      grayImage[i][j] = toByte(sum);
    }
  }

  return grayImage;
}

Float
calculateCore(int i, int j, const Matrix& image, const Matrix& core,
              PSideFunction sideFunction)
{
  int widthBorder = core[0].size() / 2;
  int heightBorder = core.size() / 2;
  Float sum = 0;

  for (int u(0); u < core.size(); ++u) {
    for (int v(0); v < core[0].size(); ++v) {
      int j1 = j - widthBorder + v;
      int i1 = i - heightBorder + u;
      if (i1 < 0 || j1 < 0 || j1 >= image[0].size() || i1 >= image.size()) {
        sum += sideFunction(image, i1, j1) * core[u][v];
      } else {
        sum += image[i1][j1] * core[u][v];
      }
    }
  }

  return sum;
}

Float
calculateWithVectors(int i, int j, const Matrix& image, const Matrix& core,
                     PSideFunction sideFunction)
{
  int widthBorder = core[0].size() / 2;
  int heightBorder = core[0].size() / 2;

  Vektor vec;

  for (int u(0); u < core[0].size(); ++u) {
    Float sum = 0;
    for (int v(0); v < core[0].size(); ++v) {
      int x = j - widthBorder + v;
      int y = i - heightBorder + u;

      if (x < 0 || y < 0 || x >= image[0].size() || y >= image.size()) {
        sum += sideFunction(image, x, y) * core[0][v];
      } else {
        sum += image[x][y] * core[0][v];
      }
    }
    vec.push_back(sum);
  }

  Float sum = 0;
  for (int u(0); u < vec.size(); ++u) {
    sum += vec[u] * core[1][u];
  }

  return sum;
}

Matrix
gaussianFilter(double sigma)
{
  Matrix res;
  int halfSize = round(3 * sigma);
  int fullSize = (halfSize << 1) + 1;
  res.fill(*new Vektor(fullSize), fullSize);
  Float r, s = 2.0 * sigma * sigma;

  Float sum = 0.0;

  double mainK = 1.0 / (M_PI * s);

  for (int x = -halfSize; x <= halfSize; x++) {
    for (int y = -halfSize; y <= halfSize; y++) {
      r = x * x + y * y;
      res[x + halfSize][y + halfSize] = mainK * expl(-(r) / s);
      sum += res[x + halfSize][y + halfSize];
    }
  }

  for (int i = 0; i < fullSize; ++i)
    for (int j = 0; j < fullSize; ++j)
      res[i][j] /= sum;

  return res;
}

Matrix
transposition(const Matrix& matrix, CoreMode coreMode)
{
  Matrix res;

  switch (coreMode) {
    case CoreMode::MONOLITH: {
      res.fill(*new Vektor(matrix.size()), matrix[0].size());
      for (int i(0); i < matrix.size(); ++i) {
        for (int j(0); j < matrix[0].size(); ++j) {
          res[j][i] = matrix[i][j];
        }
      }
    } break;
    case CoreMode::SEPARATE: {
      res = matrix;
      res[0] = matrix[1];
      res[1] = matrix[0];
    } break;
    default:
      break;
  }

  return res;
}

Matrix
normalize(const Matrix& matrix)
{
  Matrix res;
  res.fill(*new Vektor(matrix[0].size()), matrix.size());
  Vektor vecMax;
  Vektor vecMin;
  for (int i(0); i < matrix.size(); ++i) {
    vecMax.push_back(*std::max_element(matrix[i].begin(), matrix[i].end()));
    vecMin.push_back(*std::min_element(matrix[i].begin(), matrix[i].end()));
  }

  Float max = (Float)(*std::max_element(vecMax.begin(), vecMax.end()));
  Float min = (Float)(*std::min_element(vecMin.begin(), vecMin.end()));

  for (int i(0); i < matrix.size(); ++i) {
    for (int j(0); j < matrix[i].size(); ++j) {
      Float a = matrix[i][j] - min;
      Float b = max - min;
      res[i][j] = (int)((a / b) * 255);
    }
  }

  return res;
}

void
saveImage(const Matrix& matrix, const QString& fileName,
          bool normalizeFlag = true)
{
  Matrix copyMatrix = matrix;

  if (normalizeFlag) {
    copyMatrix = normalize(matrix);
  }

  QImage image(copyMatrix[0].size(), copyMatrix.size(),
               QImage::Format::Format_RGB32);

  for (int i(0); i < copyMatrix.size(); ++i) {
    for (int j(0); j < copyMatrix[0].size(); ++j) {
      int color = copyMatrix[i][j];
      image.setPixelColor(j, i, QColor(color, color, color));
    }
  }

  if (!image.save(fileName + ".jpg")) {
    outstr("File " + fileName + " not save");
  }
}

void
saveImage(const Matrix& matrix, const QString& fileName,
          const PointsVector& pVec, bool normalizeFlag = true)
{
  Matrix copyMatrix = matrix;

  if (normalizeFlag) {
    copyMatrix = normalize(matrix);
  }

  QImage image(copyMatrix[0].size(), copyMatrix.size(),
               QImage::Format::Format_RGB32);

  for (int i(0); i < copyMatrix.size(); ++i) {
    for (int j(0); j < copyMatrix[0].size(); ++j) {
      int color = copyMatrix[i][j];
      image.setPixelColor(j, i, QColor(color, color, color));
    }
  }

  for (int i(0); i < pVec.size(); ++i) {
    image.setPixelColor(pVec[i].y(), pVec[i].x(), Qt::red);
  }

  if (!image.save(fileName + ".jpg")) {
    outstr("File " + fileName + " not save");
  }
}

PSideFunction
getSideFunction(SideMode mode)
{
  PSideFunction sideFunction;
  switch (mode) {
    case SideMode::BLACK: {
      sideFunction = &blackSide;
    } break;
    case SideMode::COPY: {
      sideFunction = &boundSide;
    } break;
    case SideMode::REFLECT: {
      sideFunction = &reflectSide;
    } break;
    case SideMode::ROUND: {
      sideFunction = &roundSide;
    } break;
    default:
      break;
  }
  return sideFunction;
}

Matrix
convolution(const ImageProcessing& imgProc)
{
  Matrix res;
  res.fill(*new Vektor(imgProc.img[0].size()), imgProc.img.size());

  PSideFunction sideFunction;

  sideFunction = getSideFunction(imgProc.sideMode);

  Float (*calculateFunction)(int, int, const Matrix&, const Matrix&,
                             PSideFunction);

  switch (imgProc.coreMode) {
    case CoreMode::MONOLITH: {
      calculateFunction = &calculateCore;
    } break;
    case CoreMode::SEPARATE: {
      calculateFunction = &calculateWithVectors;
    } break;
    default:
      break;
  }

  auto coreFunction = std::bind(calculateFunction, _1, _2, imgProc.img,
                                imgProc.core, sideFunction);

  for (int i(0); i < imgProc.img.size(); ++i) {
    for (int j(0); j < imgProc.img[0].size(); ++j) {
      res[i][j] = coreFunction(i, j);
    }
  }

  return res;
}

Matrix
sobelGradient(const ImageProcessing imgProc)
{
  Matrix derX = convolution(imgProc);

  ImageProcessing imgProcInv = imgProc;
  imgProcInv.core = transposition(imgProc.core, imgProc.coreMode);

  Matrix derY = convolution(imgProcInv);

  Matrix gradient;

  int imageHeight = imgProc.img.size();
  int imageWidth = imgProc.img[0].size();

  gradient.fill(*new Vektor(imageWidth), imageHeight);

  for (int i(0); i < imageHeight; ++i) {
    for (int j(0); j < imageWidth; ++j) {
      int sum = sqrt(derX[i][j] * derX[i][j] + derY[i][j] * derY[i][j]);
      gradient[i][j] = sum;
    }
  }

  return gradient;
}

Matrix
compressImage(const ImageProcessing& imgProc)
{
  Matrix res;
  res.fill(*new Vektor(imgProc.img[0].size() >> 1), imgProc.img.size() >> 1);

  for (int i(0); i < res.size(); ++i) {
    for (int j(0); j < res[0].size(); ++j) {
      res[i][j] = imgProc.img[i << 1][j << 1];
    }
  }

  return res;
}

Matrix
gauss(const ImageProcessing imgProc)
{
  Matrix resMatrix = convolution(imgProc);

  return resMatrix;
}

bool
buildingOctave(ElementOfOctave& elementOfOctave,
               QVector<ElementOfOctave>& octave, int numberOfSpace, Float delta)
{
  double count = delta;
  for (int i(0); i < numberOfSpace - 1; ++i) {
    elementOfOctave.globalSigma *= delta;
    elementOfOctave.localSigma = count;
    elementOfOctave.imgProc.img = convolution(elementOfOctave.imgProc);

    QString currentSigma = QString::number(count);
    QString sGlobalSigma = QString::number((double)elementOfOctave.globalSigma);
    QString sNumberOfOctave = QString::number(elementOfOctave.numberOfOctave);
    QString sNumberOfImage = QString::number(i + 1);

    saveImage(elementOfOctave.imgProc.img, "octave_" + sNumberOfOctave + "_" +
                                             sNumberOfImage + "_" +
                                             currentSigma + "_" + sGlobalSigma);

    count *= delta;
    octave.push_back(elementOfOctave);
  }

  elementOfOctave.imgProc.img = compressImage(elementOfOctave.imgProc);
  return true;
}

QVector<QVector<ElementOfOctave>>
startBuildingGaussPyramid(const ImageProcessing& imgProc, int numberOfSpace)
{
  ImageProcessing tmpImgProc = imgProc;
  Float startSigma = 0.5;
  Float endSigma = 1.0;
  Float deltaS = sqrt(endSigma * endSigma - startSigma * startSigma);
  tmpImgProc.core = gaussianFilter(deltaS);
  tmpImgProc.img = convolution(tmpImgProc);

  Float globalSigma = endSigma;

  QVector<ElementOfOctave> octave;
  QVector<QVector<ElementOfOctave>> pyramid;

  ElementOfOctave elementOfOctave;
  elementOfOctave.imgProc = tmpImgProc;
  elementOfOctave.globalSigma = globalSigma;
  elementOfOctave.localSigma = 1;
  elementOfOctave.numberOfOctave = 1;

  octave.push_back(elementOfOctave);

  saveImage(tmpImgProc.img,
            "octave_1_0_1_" + QString::number((double)endSigma));

  Float mainDelta = pow(2, 1. / (numberOfSpace - 1));
  Matrix deltaCore = gaussianFilter(mainDelta);
  elementOfOctave.imgProc.core = deltaCore;

  QString sGlobalSigma;
  int currentOctave = 1;
  QString sCurrentOctave;
  bool f;

  do {
    f = buildingOctave(elementOfOctave, octave, numberOfSpace, mainDelta);
    pyramid.push_back(octave);
    octave.clear();
    ++currentOctave;

    sGlobalSigma = QString::number((double)elementOfOctave.globalSigma);
    sCurrentOctave = QString::number(currentOctave);

    if (elementOfOctave.imgProc.img.size() < 128) {
      break;
    }

    elementOfOctave.numberOfOctave = currentOctave;
    elementOfOctave.localSigma = 1;
    octave.push_back(elementOfOctave);
    saveImage(elementOfOctave.imgProc.img,
              "octave_" + sCurrentOctave + "_0_1_" + sGlobalSigma);
  } while (f);

  return pyramid;
}

double
getColorOfPointAtSigma(const QVector<QVector<ElementOfOctave>>& pyramid, int x,
                       int y, Float sigma)
{
  int numberOfSpace = pyramid[0].size();
  Float mainDelta = pow(2, 1. / (numberOfSpace - 1));
  int num = round(log(sigma) / log(mainDelta)) + 1;

  num = num + (num / numberOfSpace - 1);

  auto img = pyramid[num / (numberOfSpace)][num % (numberOfSpace)];
  int newX = x / (2 << (img.numberOfOctave - 2));
  int newY = y / (2 << (img.numberOfOctave - 2));

  std::cout << "num: " << num << std::endl;
  std::cout << "x: " << x << std::endl;
  std::cout << "y: " << y << std::endl;
  std::cout << "newX: " << newX << std::endl;
  std::cout << "newY: " << newY << std::endl;
  std::cout << "octave: " << img.numberOfOctave << std::endl;
  std::cout << "color: " << img.imgProc.img[newY][newX] << std::endl;

  return img.imgProc.img[newY][newX];
}

Float
calculateVicinity(int y, int x, const ImageProcessing& imgProc)
{
  Vektor p;
  p.fill(0, 8);
  for (int i(-1); i <= 1; ++i) {
    for (int j(-1); j <= 1; ++j) {
      Float mainC = imgProc.img[y + i][x + j];
      p[0] += (pow(imgProc.img[y + i - 1][x + j - 1] - mainC, 2));
      p[1] += (pow(imgProc.img[y + i - 1][x + j] - mainC, 2));
      p[2] += (pow(imgProc.img[y + i - 1][x + j + 1] - mainC, 2));
      p[3] += (pow(imgProc.img[y + i][x + j - 1] - mainC, 2));
      p[4] += (pow(imgProc.img[y + i][x + j + 1] - mainC, 2));
      p[5] += (pow(imgProc.img[y + i + 1][x + j - 1] - mainC, 2));
      p[6] += (pow(imgProc.img[y + i + 1][x + j] - mainC, 2));
      p[7] += (pow(imgProc.img[y + i + 1][x + j + 1] - mainC, 2));
    }
  }
  return (*std::min_element(p.begin(), p.end()));
}

PointsVector
moravec(const ImageProcessing& imgProc)
{
  Matrix mapProbability = imgProc.img;
  PointsVector res;
  for (int i(2); i < imgProc.img.size() - 2; ++i) {
    for (int j(2); j < imgProc.img[0].size() - 2; ++j) {
      mapProbability[i][j] = calculateVicinity(i, j, imgProc);
      if (mapProbability[i][j] > 1200) {
        res.push_back({ i, j });
      }
    }
  }

  return res;
}

int
main(int argc, char* argv[])
{
  QCoreApplication a(argc, argv);
  QImage mainImage;

  mainImage = loadImage("test4");

  Matrix grayImage = rgbImageToGrayscaleImage(mainImage);

  saveImage(grayImage, "gray", false);

  ImageProcessing imgProc;
  imgProc.img = grayImage;
  //  imgProc.core = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
  imgProc.core = gaussianFilter(10);
  imgProc.sideMode = SideMode::REFLECT;
  imgProc.coreMode = CoreMode::MONOLITH;

  //  Matrix sob = gauss(imgProc);
  //  saveImage(sob, "gauss_1", false);

  //  Matrix smallImg = compressImage(imgProc);
  //  saveImage(smallImg, "compress");

  //  imgProc.img = smallImg;
  //  smallImg = compressImage(imgProc);
  //  saveImage(smallImg, "compress_2");

  auto points = moravec(imgProc);
  saveImage(imgProc.img, "points", points);

  //  imgProc.img = sob;
  //  sob = gauss(imgProc);
  //  saveImage(sob, "gauss_2", false);

  //  imgProc.img = sob;
  //  sob = gauss(imgProc);
  //  saveImage(sob, "gauss_3", false);

  //  auto pyramid = startBuildingGaussPyramid(imgProc, 4);
  //  getColorOfPointAtSigma(pyramid, 16, 16, 2.5);

  outstr("All done");

  return a.exec();
}
