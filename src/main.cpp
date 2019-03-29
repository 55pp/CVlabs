#include <QCoreApplication>
#include <iostream>
#include <QBitmap>
#include <QVector>
#include <functional>
#include <algorithm>
#include <cvfunctions.h>
#include <QPainter>
#include <QSet>
#include <QMap>
#include <thread>

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

Matrix operator-(const Matrix& a, const Matrix& b)
{
  assert((a.size() == b.size()) && (a[0].size() == b[0].size()));
  Matrix c(a);
  for (int i(0); i < a.size(); ++i) {
    for (int j(0); j < a[0].size(); ++j) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }
  return c;
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

Matrix
normalizeF(const Matrix& matrix)
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
      res[i][j] = (a / b);
    }
  }

  return res;
}

Vektor
normalizeF(const Vektor& vec)
{
  Vektor res;
  res.fill(0, vec.size());

  Float max = (Float)(*std::max_element(vec.begin(), vec.end()));
  Float min = (Float)(*std::min_element(vec.begin(), vec.end()));

  for (int i(0); i < vec.size(); ++i) {
    Float a = vec[i] - min;
    Float b = max - min;
    if (b == 0) {
      continue;
    }
    res[i] = (a / b);
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

  QPainter painter(&image);
  QPen pen;
  pen.setWidth(3);
  pen.setColor(Qt::green);
  painter.setPen(pen);
  for (int i(0); i < pVec.size(); ++i) {
    painter.drawPoint(pVec[i].y(), pVec[i].x());
  }

  if (!image.save(fileName + ".jpg")) {
    outstr("File " + fileName + " not save");
  }
}

void
saveImage(const Matrix& matrix_1, const Matrix& matrix_2,
          const PointsVector& pVec_1, const PointsVector& pVec_2,
          const QVector<QPair<QPoint, QPoint>>& matchPoint,
          const QString& fileName, bool normalizeFlag = true)
{
  Matrix copyMatrix_1 = matrix_1;
  Matrix copyMatrix_2 = matrix_2;

  if (normalizeFlag) {
    copyMatrix_1 = normalize(matrix_1);
    copyMatrix_2 = normalize(matrix_2);
  }

  int imageWidth = matrix_1[0].size() + matrix_2[0].size();
  int imageHeight = std::max(matrix_1.size(), matrix_2.size());

  QImage image(imageWidth, imageHeight, QImage::Format::Format_RGB32);

  for (int i(0); i < copyMatrix_1.size(); ++i) {
    for (int j(0); j < copyMatrix_1[0].size(); ++j) {
      int color = copyMatrix_1[i][j];
      image.setPixelColor(j, i, QColor(color, color, color));
    }
  }

  int shift = copyMatrix_1[0].size();

  for (int i(0); i < copyMatrix_2.size(); ++i) {
    for (int j(0); j < copyMatrix_2[0].size(); ++j) {
      int color = copyMatrix_2[i][j];
      image.setPixelColor(j + shift, i, QColor(color, color, color));
    }
  }

  QPainter painter(&image);
  QPen pen;
  pen.setWidth(3);
  pen.setColor(Qt::green);
  painter.setPen(pen);

  for (int i(0); i < pVec_1.size(); ++i) {
    painter.drawPoint(pVec_1[i].y(), pVec_1[i].x());
  }

  for (int i(0); i < pVec_2.size(); ++i) {
    painter.drawPoint(pVec_2[i].y() + shift, pVec_2[i].x());
  }

  pen.setColor(Qt::blue);
  pen.setWidth(2);
  painter.setPen(pen);

  srand(time(NULL));

  for (int i(0); i < matchPoint.size(); ++i) {
    int x1 = matchPoint[i].first.y();
    int y1 = matchPoint[i].first.x();
    int x2 = matchPoint[i].second.y();
    int y2 = matchPoint[i].second.x();
    painter.drawLine(x1, y1, x2 + shift, y2);
    QColor color = *new QColor(qrand() % 256, qrand() % 256, qrand() % 256);
    pen.setColor(color);
    painter.setPen(pen);
  }

  if (!image.save(fileName + ".jpg")) {
    outstr("File " + fileName + " not save");
  }
}

void
saveImage(const Matrix& matrix_1, const Matrix& matrix_2,
          const PointsVector& pVec_1, const PointsVector& pVec_2,
          const QVector<QPair<QPair<QPoint, QPoint>, Float>>& matchPoint,
          const QString& fileName, bool normalizeFlag = true)
{
  Matrix copyMatrix_1 = matrix_1;
  Matrix copyMatrix_2 = matrix_2;

  if (normalizeFlag) {
    copyMatrix_1 = normalize(matrix_1);
    copyMatrix_2 = normalize(matrix_2);
  }

  int imageWidth = matrix_1[0].size() + matrix_2[0].size();
  int imageHeight = std::max(matrix_1.size(), matrix_2.size());

  QImage image(imageWidth, imageHeight, QImage::Format::Format_RGB32);

  for (int i(0); i < copyMatrix_1.size(); ++i) {
    for (int j(0); j < copyMatrix_1[0].size(); ++j) {
      int color = copyMatrix_1[i][j];
      image.setPixelColor(j, i, QColor(color, color, color));
    }
  }

  int shift = copyMatrix_1[0].size();

  for (int i(0); i < copyMatrix_2.size(); ++i) {
    for (int j(0); j < copyMatrix_2[0].size(); ++j) {
      int color = copyMatrix_2[i][j];
      image.setPixelColor(j + shift, i, QColor(color, color, color));
    }
  }

  QPainter painter(&image);
  QPen pen;
  pen.setWidth(3);
  pen.setColor(Qt::green);
  painter.setPen(pen);

  for (int i(0); i < pVec_1.size(); ++i) {
    painter.drawPoint(pVec_1[i].y(), pVec_1[i].x());
  }

  for (int i(0); i < pVec_2.size(); ++i) {
    painter.drawPoint(pVec_2[i].y() + shift, pVec_2[i].x());
  }

  pen.setColor(Qt::blue);
  pen.setWidth(2);
  painter.setPen(pen);

  srand(time(NULL));

  for (int i(0); i < matchPoint.size(); ++i) {
    int x1 = matchPoint[i].first.first.y();
    int y1 = matchPoint[i].first.first.x();
    int x2 = matchPoint[i].first.second.y();
    int y2 = matchPoint[i].first.second.x();
    painter.drawLine(x1, y1, x2 + shift, y2);
    QColor color = *new QColor(qrand() % 256, qrand() % 256, qrand() % 256);
    pen.setColor(color);
    painter.setPen(pen);
  }

  if (!image.save(fileName + ".jpg")) {
    outstr("File " + fileName + " not save");
  }
}

void
saveImage(const Matrix& matrix, const QVector<QPair<QPoint, int>>& orient,
          const QString& fileName, bool normalizeFlag = true)
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

  QPainter painter(&image);
  QPen pen;
  pen.setWidth(3);
  pen.setColor(Qt::green);
  painter.setPen(pen);

  pen.setColor(Qt::blue);
  pen.setWidth(2);
  painter.setPen(pen);

  srand(time(NULL));

  for (auto& pair : orient) {

    //    painter.rotate(pair.second);

    painter.drawPoint(pair.first.y(), pair.first.x());
    pen.setColor(Qt::blue);
    pen.setWidth(2);
    painter.setPen(pen);
    int xAfterRotate = (15) * cos(pair.second) + (0) * sin(pair.second);
    int yAfterRotate = (0) * cos(pair.second) - (15) * sin(pair.second);
    painter.drawLine(pair.first.y(), pair.first.x(),
                     pair.first.y() + xAfterRotate,
                     pair.first.x() + yAfterRotate);
    pen.setColor(Qt::green);
    pen.setWidth(3);
    painter.setPen(pen);
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
sobelDirectionOfGradient(const ImageProcessing imgProc,
                         bool inverseFlag = false)
{
  Matrix derX = convolution(imgProc);

  ImageProcessing imgProcInv = imgProc;
  imgProcInv.core = transposition(imgProc.core, imgProc.coreMode);

  Matrix derY = convolution(imgProcInv);

  Matrix direction;

  int imageHeight = imgProc.img.size();
  int imageWidth = imgProc.img[0].size();

  direction.fill(*new Vektor(imageWidth), imageHeight);

  for (int i(0); i < imageHeight; ++i) {
    for (int j(0); j < imageWidth; ++j) {
      Float sum;
      if (inverseFlag) {
        sum = atan2(derY[i][j], derX[i][j]);
      } else {
        sum = atan2(derX[i][j], derY[i][j]);
      }
      sum = sum * (180.0 / M_PI) /*+ 179*/;
      //      if (sum >= 360) {
      //        sum -= 360;
      //      }
      direction[i][j] = sum;
    }
  }

  return direction;
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

  ElementOfOctave cpyElementOfOctave = elementOfOctave;

  for (int i(0); i < 3; ++i) {
    cpyElementOfOctave.globalSigma *= delta;
    cpyElementOfOctave.localSigma = count;
    cpyElementOfOctave.imgProc.img = convolution(cpyElementOfOctave.imgProc);

    QString currentSigma = QString::number(count);
    QString sGlobalSigma =
      QString::number((double)cpyElementOfOctave.globalSigma);
    QString sNumberOfOctave =
      QString::number(cpyElementOfOctave.numberOfOctave);
    QString sNumberOfImage = QString::number(i + numberOfSpace - 1);

    saveImage(cpyElementOfOctave.imgProc.img,
              "octave_" + sNumberOfOctave + "_" + sNumberOfImage + "_" +
                currentSigma + "_" + sGlobalSigma);

    count *= delta;
    octave.push_back(cpyElementOfOctave);
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

void
checkXY(const ImageProcessing& imgProc, int y, int x, Float sourceXY,
        Vektor& vec, int index)
{
  if (x < 0 || y < 0 || x >= imgProc.img[0].size() || y >= imgProc.img.size()) {
    return;
  }

  vec[index] += (pow(imgProc.img[y][x] - sourceXY, 2));
}

PointsVector
nonMaximum(const Matrix& mapProbability, int radius, int numberOfPoints = 0)
{
  PointsVector res;

  if (numberOfPoints != 0) {
    radius = 0;
  }

  do {
    res.clear();

    for (int i(0); i < mapProbability.size(); ++i) {
      for (int j(0); j < mapProbability[0].size(); ++j) {
        Float sourceVal = mapProbability[i][j];
        if (sourceVal == 0.0) {
          continue;
        }
        bool max = true;
        for (int u(-radius); u < radius && max; ++u) {
          for (int v(-radius); v < radius && max; ++v) {
            if (!(u == 0 && v == 0)) {
              int newX = j + v;
              int newY = i + u;
              if (newX < 0 || newY < 0 || newX >= mapProbability[0].size() ||
                  newY >= mapProbability.size()) {
                continue;
              }
              Float val = mapProbability[newY][newX];
              if ((val >= sourceVal)) {
                max = false;
              }
            }
          }
        }
        if (max) {
          res.push_back({ i, j });
        }
      }
    }

    radius += 1;

  } while (res.size() > numberOfPoints && numberOfPoints != 0);

  return res;
}

Float
calculateVicinity(int y, int x, const ImageProcessing& imgProc)
{
  Vektor p;
  p.fill(0, 8);
  for (int i(-1); i <= 1; ++i) {
    for (int j(-1); j <= 1; ++j) {
      Float mainC = imgProc.img[y + i][x + j];
      checkXY(imgProc, y + i - 1, x + j - 1, mainC, p, 0);
      checkXY(imgProc, y + i - 1, x + j, mainC, p, 1);
      checkXY(imgProc, y + i - 1, x + j + 1, mainC, p, 2);
      checkXY(imgProc, y + i, x + j - 1, mainC, p, 3);
      checkXY(imgProc, y + i, x + j + 1, mainC, p, 4);
      checkXY(imgProc, y + i + 1, x + j - 1, mainC, p, 5);
      checkXY(imgProc, y + i + 1, x + j, mainC, p, 6);
      checkXY(imgProc, y + i + 1, x + j + 1, mainC, p, 7);
    }
  }
  return (*std::min_element(p.begin(), p.end()));
}

PointsVector
moravec(const ImageProcessing& imgProc, int radiusNonMax, int numberOfPoints,
        int threshold)
{
  Matrix mapProbability = imgProc.img;
  for (int i(0); i < imgProc.img.size(); ++i) {
    for (int j(0); j < imgProc.img[0].size(); ++j) {
      mapProbability[i][j] = 0;
    }
  }
  PointsVector res;
  int windowHeight = 1;
  int windowWidth = 1;

  for (int i(windowHeight); i < imgProc.img.size() - windowHeight; ++i) {
    for (int j(windowWidth); j < imgProc.img[0].size() - windowWidth; ++j) {
      mapProbability[i][j] = calculateVicinity(i, j, imgProc);
    }
  }

  saveImage(mapProbability, "map_m");
  //  mapProbability = normalize(mapProbability);

  auto vec = nonMaximum(mapProbability, radiusNonMax, numberOfPoints);
  //  res = vec;

  for (const auto& el : vec) {
    if (mapProbability[el.x()][el.y()] > threshold) {
      res.push_back(el);
    }
  }

  return res;
}

PointsVector
harris(const ImageProcessing& imgProc, int radiusOfVicinity, int radiusNonMax,
       int numberOfPoints, int threshold)
{
  ImageProcessing sobelProc = imgProc;
  sobelProc.core = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };

  Matrix derX = convolution(sobelProc);

  ImageProcessing imgProcInv = sobelProc;
  imgProcInv.core = transposition(sobelProc.core, sobelProc.coreMode);

  Matrix derY = convolution(imgProcInv);

  Matrix mapProbability = imgProc.img;
  for (int i(0); i < imgProc.img.size(); ++i) {
    for (int j(0); j < imgProc.img[0].size(); ++j) {
      mapProbability[i][j] = 0;
    }
  }

  int windowHeight = radiusOfVicinity;
  int windowWidth = radiusOfVicinity;

  Matrix gaussMatr = gaussianFilter(radiusOfVicinity / 3.);

  PointsVector res;
  for (int i(windowHeight); i < imgProc.img.size() - windowHeight; ++i) {
    for (int j(windowWidth); j < imgProc.img[0].size() - windowWidth; ++j) {
      Float a = 0;
      Float b = 0;
      Float c = 0;

      for (int u(-windowHeight); u <= windowHeight; ++u) {
        for (int v(-windowWidth); v <= windowWidth; ++v) {
          int newX = j + v;
          int newY = i + u;
          Float gaussNum =
            gaussMatr[u + radiusOfVicinity][v + radiusOfVicinity];
          a += derX[newY][newX] * derX[newY][newX] * gaussNum;
          b += derX[newY][newX] * derY[newY][newX] * gaussNum;
          c += derY[newY][newX] * derY[newY][newX] * gaussNum;
        }
      }
      //      Float det = a * c - b * b;
      //      Float trace = (a + c);
      //      Float val = det - (trace * trace * 0.05);
      //      if (trace <= 0. || val > 0) {
      //        mapProbability[i][j] = val;
      //      } else {
      //        mapProbability[i][j] = 0;
      //      }

      Float det = a * c - b * b;
      Float trace = (a + c);
      Float val = det / trace;
      if (trace > 0. || val > 0) {
        mapProbability[i][j] = val;
      }

      //      Float cc = -(b * b) + a * c;
      //      Float bb = -(a + c);
      //      Float dis = bb * bb - 4 * cc;
      //      Float l1 = (-bb + dis) / 2;
      //      Float l2 = (-bb - dis) / 2;
      //      Float lambda = std::min(l1, l2);
      //      if (lambda > threshold) {
      //        mapProbability[i][j] = lambda;
      //      } else {
      //        mapProbability[i][j] = 0;
      //      }
    }
  }

  saveImage(mapProbability, "map_h");
  //  mapProbability = normalize(mapProbability);

  auto vec = nonMaximum(mapProbability, radiusNonMax, numberOfPoints);
  //  res = vec;

  for (const auto& el : vec) {
    if (mapProbability[el.x()][el.y()] > threshold) {
      res.push_back(el);
    }
  }

  return res;
}

void
distributionOfBaskets(Float gradient, Float direction, Vektor& gisto)
{
  int numberOfBaskets = gisto.size();

  int degreeInOneBasket = 360 / numberOfBaskets;

  int mainBasket = direction / degreeInOneBasket;

  int side = ((int)direction % degreeInOneBasket) / (degreeInOneBasket >> 1);

  int sideBasket = side ? ((mainBasket + 1) % numberOfBaskets)
                        : ((mainBasket - 1) % numberOfBaskets);

  if (sideBasket == -1) {
    sideBasket = numberOfBaskets - 1;
  }

  Float procentOfGradientInMainBasket;

  int degree = direction - mainBasket * degreeInOneBasket;
  if (degree >= 0) {
    if (degree >= degreeInOneBasket / 2) {
      procentOfGradientInMainBasket =
        (degree - degreeInOneBasket / 2) / (double)(degreeInOneBasket);
      procentOfGradientInMainBasket = 1 - procentOfGradientInMainBasket;
    } else {
      procentOfGradientInMainBasket = degree / (double)(degreeInOneBasket);
      procentOfGradientInMainBasket += 0.5;
    }
  } else {
    procentOfGradientInMainBasket = degree / (double)(degreeInOneBasket);
    procentOfGradientInMainBasket -= 0.5;
  }

  gisto[mainBasket] += gradient * procentOfGradientInMainBasket;
  gisto[sideBasket] += gradient * (1. - procentOfGradientInMainBasket);
}

Matrix
fillingHistograms(const Matrix& gradient, const Matrix& direction, int x, int y,
                  int globalArea, int localArea, int numberOfBaskets,
                  double degree = 0)
{
  int numberOfHisto = (globalArea * globalArea) / (localArea * localArea);

  Matrix gaussMatr;

  if (globalArea % 2 == 0) {
    gaussMatr = gaussianFilter(((globalArea) / 6.));
  } else {
    gaussMatr = gaussianFilter(((globalArea - 1) / 6.));
  }

  Matrix histograms(numberOfHisto);
  for (auto& vek : histograms) {
    vek.fill(0, numberOfBaskets);
  }

  double cosd = cos(degree * M_PI / 180);
  double sind = sin(degree * M_PI / 180);

  for (int i(0); i < globalArea; ++i) {
    for (int j(0); j < globalArea; ++j) {

      int newX = x + j;
      int newY = y + i;

      int xAfterRotate =
        (j - globalArea / 2) * cosd + (i - globalArea / 2) * sind;
      int yAfterRotate =
        (i - globalArea / 2) * cosd - (j - globalArea / 2) * sind;

      xAfterRotate += (globalArea / 2);
      yAfterRotate += (globalArea / 2);

      if (xAfterRotate < 0 || yAfterRotate < 0 || xAfterRotate >= globalArea ||
          yAfterRotate >= globalArea) {
        continue;
      }

      int indexOfHisto = (yAfterRotate / localArea) * (globalArea / localArea) +
                         xAfterRotate / localArea;

      // outstr(QString::number(indexOfHisto));

      Float localGradient = gradient[newY][newX] * gaussMatr[i][j];
      Float localDirection =
        //        ((int)(720 + direction[newY][newX] - degree)) % 360;
        ((int)(360 + direction[newY][newX]) % 360) - degree;

      if (localDirection < 0) {
        localDirection += 360;
      } else if (localDirection >= 360) {
        localDirection -= 360;
      }

      distributionOfBaskets(localGradient, localDirection,
                            histograms[indexOfHisto]);
    }
  }

  return histograms;
}

QVector<QPair<QPoint, Matrix>>
descriptionOfPoints(const ImageProcessing& imgProc, PointsVector& points,
                    bool flagInv = false)
{
  ImageProcessing i = imgProc;
  i.core = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };

  Matrix gradient = sobelGradient(i);

  Matrix direction = sobelDirectionOfGradient(i, flagInv);

  int globalArea = 16;
  int localArea = 4;
  int numberOfBaskets = 36;

  QVector<QPair<QPoint, Matrix>> descriptions;

  auto it = std::begin(points);

  outstr("-------------------");

  QVector<QPair<QPoint, int>> orient;

  while (it != std::end(points)) {

    int x = (*it).y() - (globalArea >> 1);
    int y = (*it).x() - (globalArea >> 1);

    if (x < 0 || y < 0 || (x + globalArea >= i.img[0].size()) ||
        (y + globalArea >= i.img.size())) {
      it = points.erase(it);
      continue;
    }

    Matrix mainHistogram = fillingHistograms(
      gradient, direction, x, y, globalArea, globalArea, numberOfBaskets);

    Vektor vecOfmax;

    Float max =
      *std::max_element(mainHistogram[0].begin(), mainHistogram[0].end());

    double degreeOfRotate =
      (360 / numberOfBaskets) * mainHistogram[0].indexOf(max) +
      (180 / numberOfBaskets);

    vecOfmax.push_back(degreeOfRotate);

    orient.push_back({ (*it), degreeOfRotate });

    mainHistogram[0][mainHistogram[0].indexOf(max)] = 0;

    Float secondMax =
      *std::max_element(mainHistogram[0].begin(), mainHistogram[0].end());

    if ((secondMax / max) >= 0.8) {
      degreeOfRotate =
        (360 / numberOfBaskets) * mainHistogram[0].indexOf(secondMax);
      vecOfmax.push_back(degreeOfRotate);
    }

    outstr(QString::number((*it).y()) + " " + QString::number((*it).x()) + " " +
           QString::number(degreeOfRotate));

    for (auto& degree : vecOfmax) {
      Matrix histograms =
        fillingHistograms(gradient, direction, x, y, globalArea, localArea,
                          numberOfBaskets, degree);

      //      for (auto& vec : histograms) {
      //        vec = normalizeF(vec);
      //      }
      histograms = normalizeF(histograms);
      descriptions.push_back({ (*it), histograms });
    }

    ++it;
  }

  saveImage(i.img, orient, "orient");

  return descriptions;
}

Float
calculateLength(const Matrix& matrix_1, const Matrix& matrix_2)
{
  Float res = 0;
  for (int i(0); i < matrix_1.size(); ++i) {
    for (int j(0); j < matrix_1[i].size(); ++j) {
      res += pow((matrix_1[i][j] - matrix_2[i][j]), 2);
    }
  }

  return (sqrt(res));
}

QVector<QPair<QPair<QPoint, QPoint>, Float>>
matchingPoints(const QVector<QPair<QPoint, Matrix>>& descriptions_1,
               const QVector<QPair<QPoint, Matrix>>& descriptions_2,
               double threshold)
{
  QVector<QPair<QPair<QPoint, QPoint>, Float>> res;
  QMap<QPoint, QPair<QPair<QPoint, QPoint>, Float>> existPoints;
  QVector<QVector<QPair<QPoint, Float>>> allLengths(descriptions_1.size());

  for (int i(0); i < descriptions_1.size(); ++i) {
    QVector<QPoint> similarPoints;
    Float minLength = 999999999;
    for (int j(0); j < descriptions_2.size(); ++j) {
      Float currentLength =
        calculateLength(descriptions_1[i].second, descriptions_2[j].second);

      if (currentLength > threshold) {
        continue;
      }

      allLengths[i].push_back({ descriptions_2[j].first, currentLength });
    }
    std::sort(allLengths[i].begin(), allLengths[i].end(),
              [](const QPair<QPoint, Float>& a, const QPair<QPoint, Float>& b) {
                return a.second < b.second;
              });
  }

  bool flag = false;

  // QVector<QVector<QPair<QPoint, Float>>>
  do {
    flag = false;
    for (auto& vector : allLengths) {
      if (vector.size() == 0) {
        continue;
      }
      QPoint pnt = vector[0].first;
      Float length = vector[0].second;
      for (auto& vector_2 : allLengths) {
        if (vector_2.size() == 0) {
          continue;
        }
        if (pnt == vector_2[0].first && length > vector_2[0].second) {
          //          length = vector_2[0].second;
          vector.remove(0);
          flag = true;
          break;
        } else if (pnt == vector_2[0].first && length < vector_2[0].second) {
          vector_2.remove(0);
          flag = true;
          break;
        }
      }
    }
  } while (flag);

  res.clear();
  for (int i(0); i < allLengths.size(); ++i) {
    if (allLengths[i].size() != 0) {
      res.push_back({ { descriptions_1[i].first, allLengths[i][0].first },
                      allLengths[i][0].second });
    }
  }

  return res;
}

QVector<QPair<QPair<QPoint, QPoint>, Float>>
matchingPointsNDR(const QVector<QPair<QPoint, Matrix>>& descriptions_1,
                  const QVector<QPair<QPoint, Matrix>>& descriptions_2,
                  double threshold)
{
  QVector<QPair<QPair<QPoint, QPoint>, Float>> res;

  QVector<QVector<QPair<QPoint, Float>>> allLengths(descriptions_1.size());

  for (int i(0); i < descriptions_1.size(); ++i) {
    QPoint similarPoint = { -1, -1 };
    Float minLength = 999999999;
    for (int j(0); j < descriptions_2.size(); ++j) {
      Float currentLength =
        calculateLength(descriptions_1[i].second, descriptions_2[j].second);

      if (currentLength > threshold) {
        continue;
      }

      allLengths[i].push_back({ descriptions_2[j].first, currentLength });
    }

    std::sort(allLengths[i].begin(), allLengths[i].end(),
              [](const QPair<QPoint, Float>& a, const QPair<QPoint, Float>& b) {
                return a.second < b.second;
              });
  }

  int i = -1;
  for (auto& vec : allLengths) {
    ++i;
    if (vec.size() == 0) {
      continue;
    }
    if (vec.size() == 1) {
      res.push_back(
        { { descriptions_1[i].first, vec[0].first }, vec[0].second });
      continue;
    }
    if ((vec[0].second / vec[1].second) >= 0.8) {
      continue;
    } else {
      res.push_back(
        { { descriptions_1[i].first, vec[0].first }, vec[0].second });
      continue;
    }
  }

  //  QVector<QPair<QPoint, QPoint>> complMatch_1;
  //  QVector<QPair<QPoint, QPoint>> complMatch_2;

  //  for (int i(0); i < descriptions_1.size(); ++i) {
  //    QPoint similarPoint = { -1, -1 };
  //    Float minLength = 999999999;
  //    for (int j(0); j < descriptions_2.size(); ++j) {
  //      Float currentLength =
  //        calculateLength(descriptions_1[i].second, descriptions_2[j].second);

  //      if (currentLength > threshold) {
  //        continue;
  //      }

  //      if (minLength > currentLength) {
  //        minLength = currentLength;
  //        similarPoint = descriptions_2[j].first;
  //      }
  //    }
  //    if (similarPoint != QPoint(-1, -1)) {
  //      complMatch_1.push_back({ descriptions_1[i].first, similarPoint });
  //    }
  //  }

  //  for (int i(0); i < descriptions_2.size(); ++i) {
  //    QPoint similarPoint = { -1, -1 };
  //    Float minLength = 999999999;
  //    for (int j(0); j < descriptions_1.size(); ++j) {
  //      Float currentLength =
  //        calculateLength(descriptions_1[j].second, descriptions_2[i].second);

  //      if (currentLength > threshold) {
  //        continue;
  //      }

  //      if (minLength > currentLength) {
  //        minLength = currentLength;
  //        similarPoint = descriptions_1[j].first;
  //      }
  //    }
  //    if (similarPoint != QPoint(-1, -1)) {
  //      complMatch_2.push_back({ similarPoint, descriptions_2[i].first });
  //    }
  //  }

  //  for (const auto& pair_1 : complMatch_1) {
  //    for (const auto& pair_2 : complMatch_2) {
  //      if (pair_1 == pair_2) {
  //        res.push_back({ { pair_1.first, pair_1.second }, 0 });
  //      }
  //    }
  //  }

  return res;
}

QVector<QVector<ElementOfOctave>>
buildingDoG(const QVector<QVector<ElementOfOctave>>& pyramid)
{
  QVector<QVector<ElementOfOctave>> DoG;

  for (int i(0); i < pyramid.size(); ++i) {
    QVector<ElementOfOctave> octave;
    for (int j(0); j < pyramid[i].size() - 1; ++j) {
      ElementOfOctave elementOfOctave;
      elementOfOctave.imgProc = pyramid[i][j].imgProc;
      elementOfOctave.localSigma =
        (pyramid[i][j].localSigma + pyramid[i][j + 1].localSigma) / 2;
      elementOfOctave.globalSigma =
        (pyramid[i][j].globalSigma + pyramid[i][j + 1].globalSigma) / 2;
      elementOfOctave.imgProc.img =
        pyramid[i][j].imgProc.img - pyramid[i][j + 1].imgProc.img;
      octave.push_back(elementOfOctave);
    }
    DoG.push_back(octave);
  }

  return DoG;
}

int
main(int argc, char* argv[])
{
  QCoreApplication a(argc, argv);
  QImage mainImage;

  mainImage = loadImage("test4");

  QImage secondImage = loadImage("test4_30");

  Matrix grayImage = rgbImageToGrayscaleImage(mainImage);

  //  saveImage(grayImage, "gray", false);

  time_t start = std::time(nullptr);

  ImageProcessing imgProc;
  imgProc.img = grayImage;
  //  imgProc.core = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } };
  //  imgProc.core = gaussianFilter(0.71);
  imgProc.core = gaussianFilter(1.2);
  imgProc.sideMode = SideMode::REFLECT;
  imgProc.coreMode = CoreMode::MONOLITH;

  ImageProcessing imgProc_2;
  imgProc_2.img = rgbImageToGrayscaleImage(secondImage);
  imgProc_2.core = gaussianFilter(1.2);
  imgProc_2.sideMode = SideMode::REFLECT;
  imgProc_2.coreMode = CoreMode::MONOLITH;

  //    Matrix sob = gauss(imgProc);
  //    saveImage(sob, "gauss_1", false);

  //  Matrix smallImg = compressImage(imgProc);
  //  saveImage(smallImg, "compress");

  //  imgProc.img = smallImg;
  //  smallImg = compressImage(imgProc);
  //  saveImage(smallImg, "compress_2");

  Matrix gaussRes = gauss(imgProc);
  imgProc.img = gaussRes;

  Matrix gaussRes_2 = gauss(imgProc_2);
  imgProc_2.img = gaussRes_2;

  //  auto points = moravec(imgProc, 0, 500, 400); // test5
  //  auto points = moravec(imgProc, 1, 800);
  //  saveImage(imgProc.img, "points_m_l", points, false);

  //  Matrix gaussRes = gauss(imgProc);
  //  imgProc.img = gaussRes;
  //  auto points = harris(imgProc, 3, 3, 0);
  //  auto points_2 = harris(imgProc, 6, 4, 0, 800); // test5

  PointsVector points_2;

  auto l1 = [&imgProc, &points_2]() {
    points_2 = harris(imgProc, 7, 0, 500, 900);
    saveImage(imgProc.img, "points_h", points_2, false);
  };

  //  auto points_2 = harris(imgProc, 7, 0, 500, 800); // test4
  //  saveImage(imgProc.img, "points_h", points_2, false);

  PointsVector points_2m;

  auto l2 = [&imgProc_2, &points_2m]() {
    points_2m = harris(imgProc_2, 7, 0, 500, 900);
    saveImage(imgProc_2.img, "points_h_2m", points_2m, false);
  };

  //  auto points_2m = harris(imgProc_2, 7, 0, 500, 800);
  //  saveImage(imgProc_2.img, "points_h_2m", points_2m, false);

  std::thread th1(l1);
  std::thread th2(l2);

  th1.join();
  th2.join();

  QVector<QPair<QPoint, Matrix>> desc_1;

  auto l3 = [&desc_1, &imgProc, &points_2]() {
    desc_1 = descriptionOfPoints(imgProc, points_2, true);
  };

  //  auto desc_1 = descriptionOfPoints(imgProc, points_2);
  //  auto desc_1_1 = descriptionOfPoints(imgProc, points_2);

  QVector<QPair<QPoint, Matrix>> desc_2;

  auto l4 = [&desc_2, &imgProc_2, &points_2m]() {
    desc_2 = descriptionOfPoints(imgProc_2, points_2m, true);
  };

  std::thread th3(l3);
  std::thread th4(l4);

  //  auto desc_2 = descriptionOfPoints(imgProc_2, points_2m);
  //  auto desc_2_1 = descriptionOfPoints(imgProc_2, points_2m);

  th3.join();
  th4.join();

  //  auto mPoint_1 = matchingPoints(desc_1, desc_2, 0.9);
  auto mPoint_2 = matchingPointsNDR(desc_1, desc_2, 2000);
  //  auto mPoint_2 = matchingPoints(desc_1_1, desc_2_1, 0.8);

  saveImage(imgProc.img, imgProc_2.img, points_2, points_2m, mPoint_2,
            "test_match");

  //  auto compMatch = tempFun(mPoint_1, mPoint_2);

  //  saveImage(imgProc.img, imgProc_2.img, points_2, points_2m, compMatch,
  //            "test_match");

  //  imgProc.img = sob;
  //  sob = gauss(imgProc);
  //  saveImage(sob, "gauss_2", false);

  //  imgProc.img = sob;
  //  sob = gauss(imgProc);
  //  saveImage(sob, "gauss_3", false);

  //  auto pyramid = startBuildingGaussPyramid(imgProc, 4);
  //  getColorOfPointAtSigma(pyramid, 16, 16, 2.5);

  outstr(QString::number(std::time(nullptr) - start));

  outstr("All done");

  return a.exec();
}
