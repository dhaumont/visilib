#pragma once

/********************************************************************
 * Convex hull code from Cédric Laugerotte
 * http://student.ulb.ac.be/~claugero/3dch/
 /********************************************************************
 *   +-------------+
 *  / Description /
 * +-------------+
 *
 * This code computes the 3D convex hull of a points set.
 *
 * It is based on the code from "Computational Geometry in C"
 * (Second Edition) by Joseph O'Rourke.
 * The original code can be found on the following website
 * http://cs.smith.edu/~orourke/books/ftp.html
 *
 * This adaptation provides a class convex_hull which computes the
 * 3D convex hull of a set of vertices.
 *
 *   +-------+
 *  / Input /
 * +-------+
 *
 * double *vertices;
 * int n_vertices
 * vertices are organized as follows:
 * i-th vertex   : vertices[3*i], vertices[3*i+1], vertices[3*i+2]
 * so, the size of vertices is 3*n_vertices
 *
 *   +--------+
 *  / Output /
 * +--------+
 *
 * after calling the method "compute", the convex hull is organized
 * as follows:
 * vertices (from the class convex_hull) contains the vertices of the
 * hull. faces contains the faces defining the hull.
 * The result can be sent in a OBJ file (export_obj) or into arrays
 * for visualization.
 *
 ********************************************************************/

#pragma once

#include<vector>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

namespace ch
{
    class convex_hull_vertex;
    class convex_hull_edge;
    class convex_hull_face;

    struct vertex
    {
        vertex(double aX, double aY, double aZ) : x(aX), y(aY), z(aZ){}

        double x;
        double y;
        double z;
    };

    class convex_hull
    {
        friend class convex_hull_edge;
        friend class convex_hull_face;
    public:
        convex_hull();
        ~convex_hull();

        void compute(void);
        void addOneVertice(double x, double y, double z);
        int  get_n_vertices(void);
        int  get_n_faces(void);

        /* output */
        int get_convex_hull(double **vertices, int *n_vertices, int **faces, int *n_faces);
        void export_obj(char *filename);
        bool getVertices(std::vector<vertex>& anArray);
    private:
        void add_vertex(convex_hull_vertex *v);
        void add_edge(convex_hull_edge *e);
        void add_face(convex_hull_face *f);
        void delete_vertex(convex_hull_vertex *v);
        void delete_edge(convex_hull_edge *e);
        void delete_face(convex_hull_face *f);

        int are_collinear(convex_hull_vertex *v1, convex_hull_vertex *v2, convex_hull_vertex *v3);
        int volume_sign(convex_hull_face *f, convex_hull_vertex *v);
        int add_one(convex_hull_vertex *v);
        void clean_up(convex_hull_vertex *vnext);
        void clean_edges(void);
        void clean_faces(void);
        void clean_vertices(convex_hull_vertex *vnext);
        int double_triangle(void);
        int construct_hull(void);

        int get_vertex_index(convex_hull_vertex *v);

    private:
        convex_hull_vertex *vertices;
        convex_hull_edge *edges;
        convex_hull_face *faces;
    };

    /**********************/
    /*** convex_hull_vertex ***/
    /**********************/
    class convex_hull_vertex
    {
        friend class convex_hull;
        friend class convex_hull_face;
    public:
        convex_hull_vertex() {}
        convex_hull_vertex(double x, double y, double z);

    private:
        double          pt[3];
        convex_hull_edge   *duplicate;
        short          on_hull;
        short          processed;
        convex_hull_vertex *prev, *next;
    };

    /********************/
    /*** convex_hull_edge ***/
    /********************/
    class convex_hull_edge
    {
        friend class convex_hull;
        friend class convex_hull_face;
    public:
        convex_hull_edge(convex_hull *hull3D);

    private:
        convex_hull_face   *adj_faces[2];
        convex_hull_vertex *end_points[2];
        convex_hull_face   *new_face;
        int            to_delete;
        convex_hull_edge   *prev, *next;
    };

    /********************/
    /*** convex_hull_face ***/
    /********************/
    class convex_hull_face
    {
        friend class convex_hull;
    public:
        convex_hull_face(convex_hull *hull3D);
        convex_hull_face(convex_hull *hull3D, convex_hull_vertex *v1, convex_hull_vertex *v2, convex_hull_vertex *v3, convex_hull_face *f);
        convex_hull_face(convex_hull *hull3D, convex_hull_edge *e, convex_hull_vertex *v);

    private:
        void make_ccw(convex_hull_edge *e, convex_hull_vertex *v);

    private:
        convex_hull_edge   *edges[3];
        convex_hull_vertex *vertices[3];
        int            visible;
        convex_hull_face   *prev, *next;
    };

    inline convex_hull::convex_hull()
    {
        vertices = NULL;
        edges = NULL;
        faces = NULL;
    }

#define CH_SWAP(t,x,y) { t = x; x = y; y = t; }
#define CH_ADD( head, p )  if ( head )  { \
    p->next = head; \
    p->prev = head->prev; \
    head->prev = p; \
    p->prev->next = p; \
        																						            } \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                																																																																																																																																																																																													            else { \
            head = p; \
            head->next = head->prev = p; \
        																						            }

#define CH_DELETE( head, p ) if ( head )  { \
    if ( head == head->next ) \
    head = NULL;  \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                																																																																																																																																																																																																																																																												                else if ( p == head ) \
                head = head->next; \
                p->next->prev = p->prev;  \
                p->prev->next = p->next; \
                delete p;\
        																						            }

    inline void convex_hull::addOneVertice(double x, double y, double z)
    {
        /* init vertices */
        add_vertex(new convex_hull_vertex(x, y, z));
    }

    inline void convex_hull::add_vertex(convex_hull_vertex *v) { CH_ADD(vertices, v); }
    inline void convex_hull::add_edge(convex_hull_edge *e)   { CH_ADD(edges, e); }
    inline void convex_hull::add_face(convex_hull_face *f)   { CH_ADD(faces, f); }
    inline void convex_hull::delete_vertex(convex_hull_vertex *v) { CH_DELETE(vertices, v); }
    inline void convex_hull::delete_edge(convex_hull_edge *e)   { CH_DELETE(edges, e); }
    inline void convex_hull::delete_face(convex_hull_face *f)   { CH_DELETE(faces, f); }

    /**********/
    /* output */
    /**********/
    inline int convex_hull::get_n_vertices(void)
    {
        convex_hull_vertex *v = vertices;
        int n = 0;
        if (!vertices)
            return 0;
        else
            do { v = v->next; n++; } while (v != vertices);
        return n;
    }

    inline int convex_hull::get_n_faces(void)
    {
        convex_hull_face *f = faces;
        int n = 0;
        if (!faces)
            return 0;
        else
            do { f = f->next; n++; } while (f != faces);
        return n;
    }

    inline int convex_hull::get_vertex_index(convex_hull_vertex *v)
    {
        convex_hull_vertex *v_walk = vertices;
        int index = 0;
        do {
            if (v_walk == v)
                return index;
            index++;
            v_walk = v_walk->next;
        } while (v_walk != vertices);
        return -1;
    }

    inline int convex_hull::get_convex_hull(double **v, int *nv, int **f, int *nf)
    {
        *nv = get_n_vertices();
        *nf = get_n_faces();

        /* memory allocation */
        double *cv_vertices = (double*)malloc(3 * (*nv)*sizeof(double));
        int   *cv_faces = (int*)malloc(3 * (*nf)*sizeof(int));
        if (!cv_vertices || !cv_faces)
        {
            v = NULL; f = NULL;
            *nv = *nf = 0;
            return 0;
        }

        /* vertices */
        convex_hull_vertex *v_walk = vertices;
        int i = 0;
        do {
            cv_vertices[3 * i] = v_walk->pt[0];
            cv_vertices[3 * i + 1] = v_walk->pt[1];
            cv_vertices[3 * i + 2] = v_walk->pt[2];
            v_walk = v_walk->next;
            i++;
        } while (v_walk != vertices);

        /* faces */
        i = 0;
        convex_hull_face *f_walk = faces;
        do {
            cv_faces[3 * i] = get_vertex_index(f_walk->vertices[0]);
            cv_faces[3 * i + 1] = get_vertex_index(f_walk->vertices[1]);
            cv_faces[3 * i + 2] = get_vertex_index(f_walk->vertices[2]);
            f_walk = f_walk->next;
            i++;
        } while (f_walk != faces);

        *v = cv_vertices;
        *f = cv_faces;

        return 1;
    }

    /*****************/
    /*** Algorithm ***/
    /*****************/
    inline void convex_hull::compute(void)
    {
        double_triangle();
        construct_hull();
    }

    /* builds the initial double triangle */
    inline int convex_hull::double_triangle(void)
    {
        convex_hull_vertex *v0, *v1, *v2, *v3;
        int vol;

        /* find 3 non collinear points */
        v0 = vertices;
        while (are_collinear(v0, v0->next, v0->next->next))
        {
            if ((v0 = v0->next) == vertices)
            {
                printf("All the vertices are collinear\n");
                return 1;
            }
        }
        v1 = v0->next;
        v2 = v1->next;

        /* mark the vertices as processed */
        v0->processed = 1;
        v1->processed = 1;
        v2->processed = 1;

        /* create the two "twins" faces */
        convex_hull_face *f0, *f1 = NULL;
        f0 = new convex_hull_face(this, v0, v1, v2, f1);
        f1 = new convex_hull_face(this, v2, v1, v0, f0);

        /* link adjacent face fields */
        f0->edges[0]->adj_faces[1] = f1;
        f0->edges[1]->adj_faces[1] = f1;
        f0->edges[2]->adj_faces[1] = f1;
        f1->edges[0]->adj_faces[1] = f0;
        f1->edges[1]->adj_faces[1] = f0;
        f1->edges[2]->adj_faces[1] = f0;

        /* find a fourth, non coplanar point to form tetrahedron */
        v3 = v2->next;
        vol = volume_sign(f0, v3);
        while (!vol)
        {
            if ((v3 = v3->next) == v0)
            {
                printf("All the vertices are coplanar\n");
                return 1;
            }
            vol = volume_sign(f0, v3);
        }

        /* insure that v3 will be the first added */
        vertices = v3;

        return 0;
    }

    /*
    * construct_hull adds the vertices to the hull one at a time.
    */
    inline int convex_hull::construct_hull(void)
    {
        convex_hull_vertex *v, *vnext;

        v = vertices;
        do {
            vnext = v->next;
            if (!v->processed)
            {
                v->processed = 1;
                add_one(v);
                clean_up(vnext);
            }
            v = vnext;
        } while (v != vertices);

        return 0;
    }

    /*
    * are_collinear checks to see if the three points given are
    * collinear by checking to see if each element of the cross
    * product is zero.
    */
    inline int convex_hull::are_collinear(convex_hull_vertex *v1, convex_hull_vertex *v2, convex_hull_vertex *v3)
    {
        return
            (fabs((v3->pt[2] - v1->pt[2]) * (v2->pt[1] - v1->pt[1]) -
            (v2->pt[2] - v1->pt[2]) * (v3->pt[1] - v1->pt[1])) < 1e-9
            &&fabs((v2->pt[2] - v1->pt[2]) * (v3->pt[0] - v1->pt[0]) -
            (v2->pt[0] - v1->pt[0]) * (v3->pt[2] - v1->pt[2])) < 1e-9
            && fabs((v2->pt[0] - v1->pt[0]) * (v3->pt[1] - v1->pt[1]) -
            (v2->pt[1] - v1->pt[1]) * (v3->pt[0] - v1->pt[0])) < 1e-9);
    }

    /*
    * volume_sign returns the sign of the volume of the tetrahedron
    * determined by f and p.
    * Volume_sign is +1 iff p is on the negative side of f, where the
    * positive side is determined by the rh-rule. So the volume is
    * positive if the ccw normal to f points outside the tetrahedron.
    * The final fewer-multiplications form is due to Bob Williamson.
    */
    inline int convex_hull::volume_sign(convex_hull_face *f, convex_hull_vertex *v)
    {
        double ax, ay, az, bx, by, bz, cx, cy, cz;
        double vol;

        ax = f->vertices[0]->pt[0] - v->pt[0];
        ay = f->vertices[0]->pt[1] - v->pt[1];
        az = f->vertices[0]->pt[2] - v->pt[2];
        bx = f->vertices[1]->pt[0] - v->pt[0];
        by = f->vertices[1]->pt[1] - v->pt[1];
        bz = f->vertices[1]->pt[2] - v->pt[2];
        cx = f->vertices[2]->pt[0] - v->pt[0];
        cy = f->vertices[2]->pt[1] - v->pt[1];
        cz = f->vertices[2]->pt[2] - v->pt[2];

        vol = ax * (by*cz - bz*cy)
            + ay * (bz*cx - bx*cz)
            + az * (bx*cy - by*cx);

        if (vol > 1e-9) return  1;
        else if (vol < -1e-9) return -1;
        else                   return  0;
    }

    /*
    * add_one is passed a vertex. It first determines all faces visible
    * from that point. If none are visible then the point is marked as
    * not on hull. Next is a loop over edges. If both faces adjacent to
    * an edge are visible, then the edge is marked for deletion. If just
    * one of the adjacent faces is visible then a new face is constructed.
    */
    inline int convex_hull::add_one(convex_hull_vertex *v)
    {
        convex_hull_face *f;
        convex_hull_edge *e, *temp;
        int vol;
        int vis = 0;

        /* marks faces visible from v */
        f = faces;
        do {
            vol = volume_sign(f, v);
            if (vol < 0)
            {
                f->visible = 1;
                vis = 1;
            }
            f = f->next;
        } while (f != faces);

        /* if no faces are visible from v, then v is inside the hull */
        if (!vis)
        {
            v->on_hull = 0;
            return 0;
        }

        /* mark edges in interior of visible region for deletion.
        erect a new face based on each border edge */
        e = edges;
        do {
            temp = e->next;
            if (e->adj_faces[0] && e->adj_faces[1])
            {
                if (e->adj_faces[0]->visible && e->adj_faces[1]->visible)
                    /* e interior: mark for deletion */
                    e->to_delete = 1;
                else if (e->adj_faces[0]->visible || e->adj_faces[1]->visible)
                    /* e border: make a new face */
                    e->new_face = new convex_hull_face(this, e, v);
            }
            e = temp;
        } while (e != edges);

        return 1;
    }

    /*
    * goes through each data structure list and clears all flags
    * and NULLs out some pointers.
    */
    inline void convex_hull::clean_up(convex_hull_vertex *vnext)
    {
        clean_edges();
        clean_faces();
        clean_vertices(vnext);
    }

    /*
    * runs through the edge list and cleans up the structure.
    * If there is a newface then it will put that face in place
    * of the visible face and NULL out newface. It also deletes
    * so marked edges.
    */
    inline void convex_hull::clean_edges(void)
    {
        if (!edges)
            return;

        convex_hull_edge *e, *t;

        /* integrate the new face's into the data structure */
        /* check every edge */
        e = edges;
        do {
            if (e->new_face)
            {
                if (e->adj_faces[0]->visible) e->adj_faces[0] = e->new_face;
                else                          e->adj_faces[1] = e->new_face;
                e->new_face = NULL;
            }
            e = e->next;
        } while (e != edges);

        /* delete any edges marked for deletion */
        while (edges && edges->to_delete)
        {
            e = edges;
            delete_edge(e);
        }
        if (edges)
        {
            e = edges->next;

            do {
                if (e->to_delete)
                {
                    t = e;
                    e = e->next;
                    delete_edge(t);
                }
                else
                    e = e->next;
            } while (e != edges);
        }
    }

    /*
    * runs through the face list and deletes any face marked visible.
    */
    inline void convex_hull::clean_faces(void)
    {
        convex_hull_face *f, *t;

        while (faces && faces->visible)
        {
            f = faces;
            delete_face(f);
        }
        if (faces)
        {
            f = faces->next;
            do {
                if (f->visible)
                {
                    t = f;
                    f = f->next;
                    delete_face(t);
                }
                else
                    f = f->next;
            } while (f != faces);
        }
    }

    /*
    * runs through the vertex list and deletes the vertice
    * that are marked as processed but are not incident to
    * any undeleted edges.
    * The pointer to vnext, pvnext, is used to alter vnext
    * in construct_hull() if we are about to delete vnext.
    */
    inline void convex_hull::clean_vertices(convex_hull_vertex *vnext)
    {
        convex_hull_edge *e;
        convex_hull_vertex *v, *t;

        /* mark all vertices incident to some undeleted edge as
        on the hull */
        if (edges)
        {
            e = edges;
            do {
                e->end_points[0]->on_hull = e->end_points[1]->on_hull = 1;
                e = e->next;
            } while (e != edges);
        }
        /* delete all vertices that have been processed but are
        not on the hull */
        if (vertices)
        {
            v = vertices;
            while (vertices && vertices->processed && !vertices->on_hull)
            {
                /* if about to delete vnext, advance it first */
                if (v == vnext)
                    vnext = v->next;
                v = vertices;
                delete_vertex(v);
            }
            v = vertices->next;
            do {
                if (v->processed && !v->on_hull)
                {
                    t = v;
                    v = v->next;
                    delete_vertex(t);
                }
                else
                    v = v->next;
            } while (v != vertices);

            /* reset flags */
            v = vertices;
            do {
                v->duplicate = NULL;
                v->on_hull = 0;
                v = v->next;
            } while (v != vertices);
        }
    }

    /**********************/
    /*** convex_hull_vertex ***/
    /**********************/
    inline convex_hull_vertex::convex_hull_vertex(double x, double y, double z)
    {
        pt[0] = x;  pt[1] = y;  pt[2] = z;
        duplicate = NULL;
        on_hull = 0;
        processed = 0;
    }

    /********************/
    /*** convex_hull_edge ***/
    /********************/
    inline convex_hull_edge::convex_hull_edge(convex_hull *hull3D)
    {
        adj_faces[0] = adj_faces[1] = new_face = NULL;
        end_points[0] = end_points[1] = NULL;
        to_delete = 0;
        hull3D->add_edge(this);
    }

    /********************/
    /*** convex_hull_face ***/
    /********************/
    inline convex_hull_face::convex_hull_face(convex_hull *hull3D)
    {
        edges[0] = edges[1] = edges[2] = NULL;
        vertices[0] = vertices[1] = vertices[2] = NULL;
        visible = 0;
        hull3D->add_face(this);
    }

    inline convex_hull_face::convex_hull_face(convex_hull *hull3D, convex_hull_vertex *v1, convex_hull_vertex *v2, convex_hull_vertex *v3, convex_hull_face *f)
    {
        convex_hull_edge *e0, *e1, *e2;

        /* create edges of the initial triangle */
        if (!f)
        {
            e0 = new convex_hull_edge(hull3D);
            e1 = new convex_hull_edge(hull3D);
            e2 = new convex_hull_edge(hull3D);
        }
        else
        {
            e0 = f->edges[2];
            e1 = f->edges[1];
            e2 = f->edges[0];
        }
        e0->end_points[0] = v1;       e0->end_points[1] = v2;
        e1->end_points[0] = v2;       e1->end_points[1] = v3;
        e2->end_points[0] = v3;       e2->end_points[1] = v1;

        /* create face for triangle */
        edges[0] = e0;    edges[1] = e1;    edges[2] = e2;
        vertices[0] = v1;    vertices[1] = v2;    vertices[2] = v3;
        visible = 0;

        /* links edges to face */
        e0->adj_faces[0] = e1->adj_faces[0] = e2->adj_faces[0] = this;

        hull3D->add_face(this);
    }

    /*
    * makes a new face and two new edges between the edge and the point
    * that are passed to it.
    */
    inline convex_hull_face::convex_hull_face(convex_hull *hull3D, convex_hull_edge *e, convex_hull_vertex *v)
    {
        convex_hull_edge *new_edges[2];
        int i, j;

        /* make two new edges (if don't already exist)*/
        for (i = 0; i < 2; ++i)
        {
            /* if the edge exists, copy it into new_edges */
            new_edges[i] = e->end_points[i]->duplicate;
            if (new_edges[i]==nullptr)
            {
                /* otherwise (duplicate is NULL) */
                new_edges[i] = new convex_hull_edge(hull3D);
                new_edges[i]->end_points[0] = e->end_points[i];
                new_edges[i]->end_points[1] = v;
                e->end_points[i]->duplicate = new_edges[i];
            }
        }
        /* make the new face */
        edges[0] = e;
        edges[1] = new_edges[0];
        edges[2] = new_edges[1];
        visible = 0;
        make_ccw(e, v);

        /* set the adjacent face pointers */
        for (i = 0; i < 2; ++i)
            for (j = 0; j < 2; ++j)
                /* only the NULL link should be set to this face */
                if (!new_edges[i]->adj_faces[j])
                {
                    new_edges[i]->adj_faces[j] = this;
                    break;
                }

        hull3D->add_face(this);
    }

    /*
    * make_ccw puts the vertices in the face structure in
    * counterclockwise order.  We want to store the vertices
    * in the same order as in the visible face.  The third
    * vertex is always p.
    *
    * Although no specific ordering of the edges of a face are
    * used by the code, the following condition is maintained
    * for each face f:
    * one of the two endpoints of f->edge[i] matches f->vertex[i].
    * But note that this does not imply that f->edge[i] is between
    * f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
    */
    inline void convex_hull_face::make_ccw(convex_hull_edge *e, convex_hull_vertex *v)
    {
        convex_hull_face *fv; /* the visible face adjacent to e */
        int i;            /* index of e->end_points[0] in fv */
        convex_hull_edge *s;  /* temporary, for swapping */

        if (e->adj_faces[0]->visible) fv = e->adj_faces[0];
        else                          fv = e->adj_faces[1];

        /* set vertices[0] & vertices[1] to have the same position as do
        the corresponding  vertices of fv */
        for (i = 0; fv->vertices[i] != e->end_points[0]; ++i);

        /* orient this the same as fv */
        if (fv->vertices[(i + 1) % 3] != e->end_points[1])
        {
            vertices[0] = e->end_points[1];
            vertices[1] = e->end_points[0];
        }
        else
        {
            vertices[0] = e->end_points[0];
            vertices[1] = e->end_points[1];
            CH_SWAP(s, edges[1], edges[2]);
        }

        /* this swap is tricky. e is edges[0]. edges[1] is based on end_points[0],
        edges[2] on end_points[1]. So if e is oriented "forwards", we need to
        move  edges[1] to follow [0], because it precedes */
        vertices[2] = v;
    }

    inline  convex_hull::~convex_hull()
    {
        convex_hull_vertex *v, *vt;
        convex_hull_edge *e, *et;
        convex_hull_face *f, *ft;

        // clean vertices
        v = vertices;
        do {
            vt = v;
            v = v->next;
            delete_vertex(vt);
        } while (vertices->next != vertices);
        delete_vertex(vertices);

        // clean edges
        if (edges)
        {
            e = edges;
            do {
                et = e;
                e = e->next;
                delete_edge(et);
            } while (edges->next != edges);
            delete_edge(edges);
        }

        // clean faces
        if (faces)
        {
            f = faces;
            do {
                ft = f;
                f = f->next;
                delete_face(ft);
            } while (faces->next != faces);
            delete_face(faces);
        }
    }

    inline bool convex_hull::getVertices(std::vector<vertex>& anArray)
    {
        if (!faces)
            return false;

        int nv = get_n_vertices();

        double *cv_vertice = new double[3 * nv];

        /* vertices */
        convex_hull_vertex *v_walk = vertices;
        int i = 0;
        do {
            cv_vertice[3 * i] = v_walk->pt[0];
            cv_vertice[3 * i + 1] = v_walk->pt[1];
            cv_vertice[3 * i + 2] = v_walk->pt[2];
            v_walk = v_walk->next;
            i++;
        } while (v_walk != vertices);

        /* faces */
        convex_hull_face *f_walk = faces;
        do {
            int myIndice1 = get_vertex_index(f_walk->vertices[0]);
            int myIndice2 = get_vertex_index(f_walk->vertices[1]);
            int myIndice3 = get_vertex_index(f_walk->vertices[2]);

            anArray.push_back(vertex(cv_vertice[myIndice1 * 3], cv_vertice[myIndice1 * 3 + 1], cv_vertice[myIndice1 * 3 + 2]));
            anArray.push_back(vertex(cv_vertice[myIndice2 * 3], cv_vertice[myIndice2 * 3 + 1], cv_vertice[myIndice2 * 3 + 2]));
            anArray.push_back(vertex(cv_vertice[myIndice3 * 3], cv_vertice[myIndice3 * 3 + 1], cv_vertice[myIndice3 * 3 + 2]));

            f_walk = f_walk->next;
        } while (f_walk != faces);

        delete [] cv_vertice;

        return true;
    }

#undef CH_SWAP
#undef CH_ADD
#undef CH_DELETE
}