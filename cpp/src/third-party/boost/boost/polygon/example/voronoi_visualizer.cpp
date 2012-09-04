// Boost.Polygon library voronoi_visualizer.cpp file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#include <vector>

#include <QtOpenGL/QGLWidget>
#include <QtGui/QtGui>

#include <boost/polygon/voronoi_builder.hpp>
#include <boost/polygon/voronoi_diagram.hpp>
#include <boost/polygon/voronoi_utils.hpp>
using namespace boost::polygon;

class GLWidget : public QGLWidget {
  Q_OBJECT
public:
  GLWidget(QMainWindow *parent = NULL) :
      QGLWidget(QGLFormat(QGL::SampleBuffers), parent),
      primary_edges_only_(false),
      internal_edges_only_(false) {
    startTimer(40);
  }

  QSize sizeHint() const {
    return QSize(600, 600);
  }

  void build(QString file_path) {
    brect_.clear();
    vb_.clear();
    vd_.clear();

    // Open file.
    QFile data(file_path);
    if (!data.open(QFile::ReadOnly)) {
      QMessageBox::warning(this, tr("Voronoi Visualizer"),
                           tr("Disable to open file ") + file_path);
    }

    // Read points from the file.
    QTextStream in_stream(&data);
    int num_point_sites = 0;
    int num_edge_sites = 0;
    int x1, y1, x2, y2;
    in_stream >> num_point_sites;
    for (int i = 0; i < num_point_sites; ++i) {
      in_stream >> x1 >> y1;
      vb_.insert_point(x1, y1);
      brect_.update(x1, y1);
    }
    in_stream >> num_edge_sites;
    for (int i = 0; i < num_edge_sites; ++i) {
      in_stream >> x1 >> y1 >> x2 >> y2;
      vb_.insert_segment(x1, y1, x2, y2);
      brect_.update(x1, y1);
      brect_.update(x2, y2);
    }
    in_stream.flush();

    // Build voronoi diagram.
    vb_.construct(&vd_);
    brect_ = voronoi_utils<coordinate_type>::scale(brect_, 1.4);
    shift_ = point_type((brect_.x_min() + brect_.x_max()) * 0.5,
                        (brect_.y_min() + brect_.y_max()) * 0.5);

    for (const_edge_iterator it = vd_.edges().begin();
        it != vd_.edges().end(); ++it) {
      if (!it->is_finite()) {
        remove_exterior(&(*it));
      }
    }

    // Update view.
    update_view_port();
  }

  void show_primary_edges_only() {
    primary_edges_only_ ^= true;
  }

  void show_internal_edges_only() {
    internal_edges_only_ ^= true;
  }

protected:
  void initializeGL() {
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
  }

  void paintGL() {
    qglClearColor(QColor::fromRgb(255, 255, 255));
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Draw voronoi sites.
    {
      glColor3f(0.0f, 0.5f, 1.0f);
      glPointSize(9);
      glBegin(GL_POINTS);
      for (const_cell_iterator it = vd_.cells().begin();
          it != vd_.cells().end(); ++it) {
        if (!it->contains_segment()) {
          glVertex2f(it->point0().x() - shift_.x(),
                     it->point0().y() - shift_.y());
        }
      }
      glEnd();
      glPointSize(6);
      glLineWidth(2.7f);
      glBegin(GL_LINES);
      for (const_cell_iterator it = vd_.cells().begin();
          it != vd_.cells().end(); ++it) {
        if (it->contains_segment()) {
          glVertex2f(it->point0().x() - shift_.x(),
                     it->point0().y() - shift_.y());
          glVertex2f(it->point1().x() - shift_.x(),
                     it->point1().y() - shift_.y());
        }
      }
      glEnd();
      glLineWidth(1.0);
    }

    // Draw voronoi vertices.
    /*{
      glColor3f(0.0f, 0.0f, 0.0f);
      glBegin(GL_POINTS);
      for (const_vertex_iterator it = vd_.vertices().begin();
          it != vd_.vertices().end(); ++it) {
        glVertex2f(it->vertex().x() - shift_.x(),
                   it->vertex().y() - shift_.y());
      }
      glEnd();
    }*/

    // Draw voronoi edges.
    {
      glColor3f(0.0f, 0.0f, 0.0f);
      glLineWidth(1.7f);
      glBegin(GL_LINES);
      for (const_edge_iterator it = vd_.edges().begin();
          it != vd_.edges().end(); ++it) {
        if (primary_edges_only_ && !it->is_primary()) {
          continue;
        }
        if (internal_edges_only_ && it->color()) {
          continue;
        }
        std::vector<point_type> vec;
        if (!it->is_finite())
          voronoi_utils<coordinate_type>::clip(*it, brect_, vec);
        else {
          coordinate_type max_error = 1E-3 * (brect_.x_max() - brect_.x_min());
          voronoi_utils<coordinate_type>::discretize(*it, max_error, vec);
        }
        for (int i = 0; i < static_cast<int>(vec.size()) - 1; ++i) {
          glVertex2f(vec[i].x() - shift_.x(), vec[i].y() - shift_.y());
          glVertex2f(vec[i+1].x() - shift_.x(), vec[i+1].y() - shift_.y());
        }
      }
      glEnd();
    }
  }

  void resizeGL(int width, int height) {
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);
  }

  void timerEvent(QTimerEvent *) {
    update();
  }

private:
  typedef double coordinate_type;
  typedef detail::point_2d<double> point_type;
  typedef voronoi_diagram<double> VD;
  typedef VD::edge_type edge_type;
  typedef VD::cell_container_type cell_container_type;
  typedef VD::cell_container_type vertex_container_type;
  typedef VD::edge_container_type edge_container_type;
  typedef VD::const_cell_iterator const_cell_iterator;
  typedef VD::const_vertex_iterator const_vertex_iterator;
  typedef VD::const_edge_iterator const_edge_iterator;

  static const std::size_t VISITED = 1;

  void remove_exterior(const VD::edge_type* edge) {
    if (edge->color())
      return;
    edge->color(VISITED);
    edge->twin()->color(VISITED);
    const voronoi_diagram<double>::vertex_type* v = edge->vertex1();
    if (v == NULL || !edge->is_primary())
      return;
    const voronoi_diagram<double>::edge_type* e = v->incident_edge();
    do {
      remove_exterior(e);
      e = e->rot_next();
    } while (e != v->incident_edge());
  }

  void update_view_port() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(brect_.x_min() - shift_.x(), brect_.x_max() - shift_.x(),
            brect_.y_min() - shift_.y(), brect_.y_max() - shift_.y(),
            -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
  }

  bounding_rectangle<double> brect_;
  point_type shift_;
  default_voronoi_builder vb_;
  voronoi_diagram<coordinate_type> vd_;
  bool primary_edges_only_;
  bool internal_edges_only_;
};

class MainWindow : public QWidget {
  Q_OBJECT
public:
  MainWindow() {
    glWidget_ = new GLWidget();
    file_dir_ = QDir(QDir::currentPath(), tr("*.txt"));
    file_name_ = tr("");

    QHBoxLayout *centralLayout = new QHBoxLayout;
    centralLayout->addWidget(glWidget_);
    centralLayout->addLayout(create_file_layout());
    setLayout(centralLayout);

    update_file_list();
    setWindowTitle(tr("Voronoi Visualizer"));
    layout()->setSizeConstraint(QLayout::SetFixedSize);
  }

private slots:
  void primary_edges_only() {
    glWidget_->show_primary_edges_only();
  }

  void internal_edges_only() {
    glWidget_->show_internal_edges_only();
  }

  void browse() {
    QString new_path = QFileDialog::getExistingDirectory(
        0, tr("Choose Directory"), file_dir_.absolutePath());
    if (new_path.isEmpty()) {
      return;
    }
    file_dir_.setPath(new_path);
    update_file_list();
  }

  void build() {
    file_name_ = file_list_->currentItem()->text();
    QString file_path = file_dir_.filePath(file_name_);
    message_label_->setText("Building...");
    glWidget_->build(file_path);
    message_label_->setText("Double click the item to build voronoi diagram:");
    setWindowTitle(tr("Voronoi Visualizer - ") + file_path);
  }

  void print_scr() {
    if (!file_name_.isEmpty()) {
      QImage screenshot = glWidget_->grabFrameBuffer(true);
      QString output_file = file_dir_.absolutePath() + tr("/") +
          file_name_.left(file_name_.indexOf('.')) + tr(".png");
      screenshot.save(output_file, 0, -1);
    }
  }

private:
  QGridLayout *create_file_layout() {
    QGridLayout *file_layout = new QGridLayout;

    message_label_ = new QLabel("Double click item to build voronoi diagram:");

    file_list_ = new QListWidget();
    file_list_->connect(file_list_,
                        SIGNAL(itemDoubleClicked(QListWidgetItem*)),
                        this,
                        SLOT(build()));

    QCheckBox *primary_checkbox = new QCheckBox("Show primary edges only.");
    connect(primary_checkbox, SIGNAL(clicked()),
        this, SLOT(primary_edges_only()));

    QCheckBox *internal_checkbox = new QCheckBox("Show internal edges only.");
    connect(internal_checkbox, SIGNAL(clicked()),
        this, SLOT(internal_edges_only()));

    QPushButton *browse_button =
        new QPushButton(tr("Browse Input Directory"));
    connect(browse_button, SIGNAL(clicked()), this, SLOT(browse()));
    browse_button->setMinimumHeight(50);

    QPushButton *print_scr_button = new QPushButton(tr("Make Screenshot"));
    connect(print_scr_button, SIGNAL(clicked()), this, SLOT(print_scr()));
    print_scr_button->setMinimumHeight(50);

    file_layout->addWidget(message_label_, 0, 0);
    file_layout->addWidget(file_list_, 1, 0);
    file_layout->addWidget(primary_checkbox, 2, 0);
    file_layout->addWidget(internal_checkbox, 3, 0);
    file_layout->addWidget(browse_button, 4, 0);
    file_layout->addWidget(print_scr_button, 5, 0);

    return file_layout;
  }

  void update_file_list() {
    QFileInfoList list = file_dir_.entryInfoList();
    file_list_->clear();
    if (file_dir_.count() == 0) {
      return;
    }
    QFileInfoList::const_iterator it;
    for (it = list.begin(); it != list.end(); it++) {
      file_list_->addItem(it->fileName());
    }
    file_list_->setCurrentRow(0);
  }

  QDir file_dir_;
  QString file_name_;
  GLWidget *glWidget_;
  QListWidget *file_list_;
  QLabel *message_label_;
};

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  MainWindow window;
  window.show();
  return app.exec();
}

#include "voronoi_visualizer.moc"
