using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace SurfNet
{
    using static DebugLog;
    using static Mathex;

    public class FaceToTriangleIdxMap : Dictionary<int, int>
    {
    }

    public partial class KineticTriangulation
    {
        public void initialize(BasicInput input, WavefrontEdgeList p_wavefront_edges, int p_restrict_component)
        {
            //   DBG_FUNC_BEGIN(DBG_KT_SETUP);

            assert(!initialized);
            wavefront_edges = p_wavefront_edges;
            restrict_component_ = p_restrict_component;
           
            skeleton.vertices.Capacity= input.Vertices.Count;

            foreach (var bv in input.Vertices)
            {
                var sv = skeleton.new_vertex(bv.p, 0.0);
            }




            //DBG(DBG_KT_SETUP) << "Input has " << ct.max_component() + 1 << " component(s): 0 .. " << ct.max_component();
            //DBG(DBG_KT_SETUP) << "Requested straight skeleton of component " << restrict_component_;

            FaceToTriangleIdxMap face_to_triangle_idx = new FaceToTriangleIdxMap();

            TriangleOriginalVertexIndexList triangle_original_vertex_indices = new TriangleOriginalVertexIndexList();

            /* initialze the basic triangulation data structure, setting up neighborhoods */
            //  int num_initial_triangles = get_num_initial_triangles(ct);
            // DBG(DBG_KT_SETUP) << "  Have " << num_initial_triangles << " initial kinetic triangles.";
            initialize_tds(input, face_to_triangle_idx, triangle_original_vertex_indices);
            create_supporting_lines(input, face_to_triangle_idx, triangle_original_vertex_indices);
            create_kinetic_vertices(input, triangle_original_vertex_indices);
            /* until here, triangle_original_vertex_indices is consistent.
             * create_bevels is the first that may flip things around.
             */
            create_bevels(input, triangle_original_vertex_indices);


            

            foreach (var v in vertices)
            {
                var sh = skeleton.new_edge();
                var next_v = vertices[ContourVertices[v.Id].NextInLAV.ID];
                var he_v = v.Halfedge;
                var he_nv = next_v.Halfedge;
                sh.Vertex =v.Vertex;
                he_v.Opposite.Next = sh;
                sh.Next = he_nv;
            }

            foreach (var v in vertices)
            {
                var contourHe = v.Halfedge.Opposite.Next;
                var f = skeleton.new_face(v.Halfedge.Opposite.Next);
                contourHe.Prev.Face = 
                contourHe.Next.Face = f;
                if (v.Halfedge.Opposite.Next.Vertex != v.Vertex)
                {

                }
            }



            store_initial_wavefront_vertices();

            tidx_in_check_refinement = new bool[triangles.Count];
            refine_triangulation_initial();

            assert_valid(-1, CORE_ZERO);
            initialized = true;

            // DBG_FUNC_END(DBG_KT_SETUP);
        }

        private void initialize_tds(BasicInput input, Dictionary<int, int> face_to_triangle_idx, TriangleOriginalVertexIndexList triangle_original_vertex_indices)
        {
            //DBG_FUNC_BEGIN(//DBG_KT_SETUP);

            // build kinetic triangulation
            int num_t = input.Triangles.Count;
            // at each input vertex we start one kinetic vertex per sector, and the number of sectors equals the vertex degree.
            // Additionally, we'll split vertices for beveling (including degree-1 vertices)

            int num_events = num_t;
            var num_v = input.get_total_degree() + 1 + num_events + input.get_num_extra_beveling_vertices();

            /** Allocate sufficient space for our arrays */
            this.vertices = new WavefrontVertexList();
            this.vertices.Capacity = (num_v + num_t);
           

            triangles.Capacity = (num_t);

            ///   vertices.Add(WavefrontVertex.make_infinite_vertex());
            WavefrontVertex infinite = null; /// vertices.back();

            /** Create kinetic triangles.  Without any vertices initially.
             * The only exception is that we already point to the correct infinite vertex.
             * Also set up the CT-face to KT index map.
             */
            //DBG(//DBG_KT_SETUP) << "#faces: " << num_t;

            foreach (var fit in input.Triangles)
            {
                assert(triangles.Count == fit.Id);
                var t = new KineticTriangle(triangles.Count,/* fit.ID*/ 0);
                assert(fit.Id == t.Id);
                triangles.Add(t);

                int i = 0;
                foreach (var he in fit.Vertices)
                {
                    triangle_original_vertex_indices.Add(he.Id);
                }
                assert(triangles.size() * 3 == triangle_original_vertex_indices.size());
                face_to_triangle_idx.Add(fit.Id, t.Id);
                //DBG(//DBG_KT_SETUP) << "Added " << &triangles.back()
                //  << "; input vidx: "
                //  << fit.vertex(0).info().original_vertex_idx << ","
                //  << fit.vertex(1).info().original_vertex_idx << ","
                //  << fit.vertex(2).info().original_vertex_idx;
            }
            /** Set up neighborhood references between triangles that are not
             * seperated by a constraint.
             */
            foreach (var fit in input.Triangles)
            {
                // if (!fit.info().matches_component(restrict_component_)) continue;

                KineticTriangle?[] n = fit.Neighbors.Select(e => e == null ? null : triangles[face_to_triangle_idx[e.Id]]).ToArray();

                var t = triangles[fit.Id];
                //DBG(//DBG_KT_SETUP) << "set neighbors for " << t << " to "
                //<< n[0] << ", "
                //<< n[1] << ", "
                //<< n[2] << ".";
                t.set_neighbors(n);
            }

            //   assert(triangles.size() == num_initial_triangles);
            //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        internal void create_supporting_lines(
            BasicInput input,
             FaceToTriangleIdxMap face_to_triangle_idx,
                TriangleOriginalVertexIndexList triangle_original_vertex_indices
            )
        {
            // DBG_FUNC_BEGIN(DBG_KT_SETUP);
            /* total number of wavefront edges during the entire propragation.
             * this is an upper bound.
             *
             * num_t is a really rough upper bound on the number of split events.
             */
            int num_t = input.Triangles.Count + input.get_num_extra_beveling_vertices();
            var num_wavefront_edges = input.get_total_degree() * 2 + num_t * 2;
            this.wavefront_edges = new WavefrontEdgeList();
            //  this.wavefront_edges.Capacity = num_wavefront_edges;

            foreach (var fit in input.Triangles)
            {
                //if (!fit.info().matches_component(restrict_component_)) continue;
                WavefrontEdge3 w = new WavefrontEdge3();

                var triangle_idx = face_to_triangle_idx[fit.Id];
                var t = triangles[triangle_idx];

                for (int i = 0; i < 3; ++i)
                {
                    if (fit.Neighbors[i] == null)
                    {
                        int vidxu = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, ccw(i));
                        int vidxv = get_basic_vertex_idx_from_triangle_vertex_indices(input, triangle_original_vertex_indices, triangle_idx, cw(i));
                        BasicVertex u = input.Vertices[vidxu];
                        BasicVertex v = input.Vertices[vidxv];

                        if (!input.has_edge(vidxu, vidxv))
                        {
                            throw new Exception(
                                $"Cannot find edge for triangulation constraint ({vidxu}, {vidxv}) in input.  Input is not a PSLG." +
                                $" v{vidxu}:{input.Vertices[vidxu].p}\n v{vidxv}:{input.Vertices[vidxv].p}" +
                                " Probably one of these points is on the interior of another edge."
                                );
                        }

                        ////   DBG(DBG_KT_SETUP) << "at " << t << "/" << i << ": vertices are " << vidxu << " and " << vidxv;

                        //WavefrontEdge buddy_wavefront_edge = null;

                        ////  var edge = input.get_edge(vidxu, vidxv);

                        //// DBG(DBG_KT_SETUP) << "at " << t << "/" << i << ": vertices are " << vidxu << " and " << vidxv;

                        //var ct_n_hdl = fit.Neighbors[i.ccw()];
                        //if (ct_n_hdl!=null &&!face_to_triangle_idx.ContainsKey(ct_n_hdl.Id))
                        //{
                        //    /* neighbor not in map */
                        //    //assert(restrict_component_ >= 0);
                        //    //assert(!ct_n_hdl.info().matches_component(restrict_component_));
                        //}
                        //else
                        //{
                        //    KineticTriangle n = triangles[face_to_triangle_idx[ct_n_hdl.Id]];
                        //    //  assert(ct_n_hdl.info().matches_component(restrict_component_));
                        //    int idx_in_n = ct_n_hdl.index(fit);
                        //    buddy_wavefront_edge = n.wavefront(idx_in_n); /* may still be null, and that's OK */
                        //}

                       // SkeletonDCELFace skeleton_dcel_face = skeleton.setup_new_input_edge(null /*buddy_wavefront_edge*/ /* may be NULL and that's OK */);

                        wavefront_edges.Add(new WavefrontEdge(u.p, v.p, /* edge.weight*/ 1, t, null));
                        w[i] = wavefront_edges.Last();
                    }
                    else
                    {
                        w[i] = null;
                    }
                }

                t.set_wavefronts(w[0], w[1], w[2]);
            }
            // DBG_FUNC_END(DBG_KT_SETUP);
        }

        private void create_kinetic_vertices(BasicInput input, TriangleOriginalVertexIndexList triangle_original_vertex_indices)
        {

            
            this.vertices.AddRange(input.Vertices.Select(x =>(WavefrontVertex) null));

            ////DBG_FUNC_BEGIN(//DBG_KT_SETUP);

            /* We iterate over all triangles, and in each triangle over each
             * vertex v.  If the triangle is the start of a clockwise triangle
             * fan about v, that is, if the correct edge incident to v is contraint,
             * we set up the kinetic vertex.
             */
            foreach (var t_it in triangles)
            {
                /* Invariant: all previous triangles have vertices set up on
                 * the "head" side of each constraint they have.
                 *
                 * They not necessarily have vertices set up on the "tail" side of
                 * each constraint or at vertices with no incident (in this triangle)
                 * constraint.
                 */
                int t_idx = t_it.Id;
                //assert(t_idx == t_it - triangles.begin());
                //DBG(//DBG_KT_SETUP) << "setting up kinetic vertices for " << &*t_it;
                for (int i = 0; i < 3; ++i)
                {
                    //DBG(//DBG_KT_SETUP) << "vertex at idx " << i;
                    if (!t_it.is_constrained(ccw(i)))
                    {
                        continue;
                    }

                    /* the edge cw of i (opposite i's ccw vertex) is constrained. */
                    AroundVertexIterator faces_it = incident_faces_iterator(t_it, i);
                    assert(faces_it.next_triangle_ccw() == null);
                    AroundVertexIterator most_cw_triangle = faces_it.most_cw();

                    WavefrontEdge l = t_it.wavefront(ccw(i));
                    WavefrontEdge r = most_cw_triangle.t().wavefront(cw(most_cw_triangle.v_in_t_idx()));

                    BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t_idx, i);

                    if (bv.degree == 1)
                    {
                        continue; // Need to split/bevel in any case.
                    }

                    assert(bv.degree < 1);

                    // XXX do not make a vertex if this is reflex and the reflex_beveling_add is greater than 1.
                    WavefrontVertex v = vertices.make_initial_vertex(skeleton, bv.Id, bv.p, l, r, false);



                    // Should not have created a vertex that needed beveling.
                    assert(bv.reflex_beveling_add == 0 || v.is_convex_or_straight());

                    for (AroundVertexIterator it = incident_faces_iterator(t_it, i); !it.IsEnd; it = it.walk_dir_cw)
                    {
                        it.t().set_vertex(it.v_in_t_idx(), v);
                        invalidate_basic_vertex_idx_in_triangle_vertex_indices(input, triangle_original_vertex_indices, it.t().Id, it.v_in_t_idx());
                     //   Log($"  setting vertex to {v.Id} in T{it.t().Id} ({it.v_in_t_idx()})");
                    };

                    //Log($"  setting vertex 1 to {v} in {l.Id}");
                    //Log($"  setting vertex 0 to {v} in {r.Id}");
                    assert(l.vertex(1) == v);
                    assert(r.vertex(0) == v);
                }
            }
            foreach (var kv in vertices)
            {
                var cv = ContourVertices[kv.Id];
                kv.NextInLAV = vertices[cv.NextInLAV.ID];
                kv.PrevInLAV = vertices[cv.PrevInLAV.ID];
            }
            //DBG_FUNC_END(//DBG_KT_SETUP);
        }

        /** mark the original vertex as unusable.
        *
        * once we have created a kinetic vertex for a basic (input) vertex,
        * we no longer want to touch the basic vertex from input.vertices().
        *
        * This function marks the entry in the triangle-vertex to original-vertex map
        * as invalid, preventing further accesses.
        */

        private static partial void invalidate_basic_vertex_idx_in_triangle_vertex_indices(
          BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          int t_idx,
          int i
        )
        {
            int idx = t_idx * 3 + i;
            assert(idx < triangle_original_vertex_indices.size());
            assert(triangle_original_vertex_indices[idx] >= 0);
            triangle_original_vertex_indices[idx] = -1;
        }

        /** Create bevels at one vertex.
        *
        * .) find the neighboring constraints,
        * .) create beveling wavefront supporting lines
        * .) create the kinetic vertices between adjacent wavefronts, and
        * .) split one or more triangles and set the vertices accordingly.
        */

        private void create_bevels_at_vertex(BasicInput input,
        TriangleOriginalVertexIndexList triangle_original_vertex_indices, KineticTriangle t, int i)

        {
            throw new NotImplementedException();
           
            /* We should not have any triangles with undefined vertices after
             * we went past the original size.  In particular, since the
             * triangle_vertex_handle_idx array is not defined for the new ones,
             * and in fact may have become invalid for some of the existing entries
             * too.
             */
           // DBG(DBG_KT_SETUP) << t << "; vertex " << i;
            WavefrontEdge l, r;

            AroundVertexIterator most_ccw_triangle = incident_faces_iterator(t, i).most_ccw();
            AroundVertexIterator most_cw_triangle = incident_faces_iterator(t, i).most_cw();
            l = most_ccw_triangle.t().wavefront(ccw(most_ccw_triangle.v_in_t_idx()));
            r = most_cw_triangle.t().wavefront(cw(most_cw_triangle.v_in_t_idx()));

         //   DBG(DBG_KT_SETUP) << " incident edges: " << *l << ", " << *r;

            /* Build the list of wavefront edges */
            /*************************************/
            List<WavefrontEdge> edges = new List<WavefrontEdge>();
            edges.push_back(l);
            BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t.Id, i);

            if (bv.degree == 1)
            {
                if (bv.reflex_beveling_add >= 2)
                {// Beveling reflex vertex

                    throw new Exception ( "Beveling degree-one vertex with 2 or more extra vertices not implemented yet" );

                }
                {
                    Log($"l:{l}");
                    Log($"r:{r}");
                    assert(l.l().l == r.l().l.Opposite());

                  if (! Mathex.AreNear(l.l().weight, r.l().weight))
                    {
                        throw new Exception ( "Unclear what to do when the incidents weights do not match in beveling" );
                          // Beveling reflex vertex
                        
                    }

                    var w = new WavefrontEdge(new WavefrontSupportingLine(l.l().l.Opposite().Perpendicular(bv.p), l.l().weight));
                    Log($"New wavefront edgew:{w}");
                    wavefront_edges.Add(w);
                    edges.push_back(w);
                    // DBG(DBG_KT_SETUP) << " New wavefront edge: " << *w;
                }
            }
            else
            {
                throw new NotImplementedException("Beveling reflex vertices not implemented yet");
            }
            edges.Add(r);

            List<WavefrontVertex> new_wavefront_vertices = new List<WavefrontVertex>();
            assert(edges.size() > 2);

           // DBG(DBG_KT_SETUP) << " Edges:";
            foreach (var e in  edges)
            {
              //  DBG(DBG_KT_SETUP) << "  - " << *e;
            }

            /* Build the list of wavefront vertices */
            /***************************************/

            
           
            var previous_edge = edges.First();
            WavefrontVertex prev_vertex = null;
            //DBG(DBG_KT_SETUP) << " New wavefront vertices:";
            foreach (var  edge in edges.Skip(1))
            {
                WavefrontVertex v = vertices.make_initial_vertex(skeleton, bv.Id,bv.p, previous_edge, edge, true);
                Log($"add {v.Debug()}");
                if (prev_vertex!=null)
                {
                    v.link_tail_to_tail(prev_vertex);
                }
                prev_vertex = v;

                new_wavefront_vertices.Add(v);
                
                previous_edge = edge;
            }

           
            /* And update the triangulation, splitting vertices */
            /****************************************************/
            var new_v_idx =0; 
            var new_v = new_wavefront_vertices.First();
            var current_edge_idx = 1;
            var current_edge = edges[current_edge_idx];
            
            var t_it = most_ccw_triangle.Copy();
            while (!t_it.IsEnd)
            {
                int t_idx = t_it.t().Id;
          
                // if we split in the iteration before, this is a duplicate (as we set it in split_vertex()), but who cares :) */
                t_it.t().set_vertex(t_it.v_in_t_idx(), new_v);
                Log($" Setting vertex to {new_v} in t:{t_it.t().Id}({t_it.v_in_t_idx()}");

                // If the normal of the right edge points into the current triangle,
                // we move to the next wavefront/vertex.  We split the triangle by
                // duplicating the triangulation vertex here, assigning the previous
                // and new wavefront vertex to each of the triangulation vertices from
                // the split.
                
                bool split_this;
                bool need_flip = false;
                assert(current_edge == (new_v).incident_wavefront_edge(1));

                if (new_v_idx == new_wavefront_vertices.Count-1)
                {
                    Log($"   At last vertex, nothing to split anymore");
                    split_this = false;
                }
                else if (t_it.next_triangle_cw() == null)
                {
                    Log($"  At last triangle, we need to split here.");
                    split_this = true;
                }
                else
                {
                    Log($" Checking if we need to split in t{t_it.t().Id}");
                    Point2 pos_t_ccw;
                    bool infinite_ccw_vertex;
                    (pos_t_ccw, infinite_ccw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, t_it.t(), ccw(t_it.v_in_t_idx()));

                    if (!infinite_ccw_vertex)
                    {
                        Log($"  All nice and finite");
                        Point2 pos_plus_normal = bv.p + (current_edge).l().normal_direction.perpendicular(OrientationEnum.RIGHT_TURN);

                        Log($"   p     {bv.p}");
                        Log($"   p+    {pos_plus_normal}");
                        Log($"   ptccw {pos_t_ccw}");
                        var orientation = Mathex.orientation(bv.p, pos_plus_normal, pos_t_ccw);
                        // DBG(DBG_KT_SETUP) << "   orientation " << orientation;
                        split_this = (orientation == OrientationEnum.LEFT_TURN);
                    }
                    else
                    {
                        Log($"  Unbounded triangle: ccw vertex is the infinite vertex");
                        //WavefrontVertex const * const cw_v = t_it.t().vertex(cw(t_it.v_in_t_idx()));
                        Log($"   cw_vertex point              {t_it.t().vertex(cw(t_it.v_in_t_idx()))}");

                        var next_t = t_it.Copy();
                        next_t.next();
                        Point2 next_pos_t_ccw;
                        bool next_infinite_ccw_vertex;
                        (next_pos_t_ccw, next_infinite_ccw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, next_t.t(), ccw(next_t.v_in_t_idx()));
                        assert(!next_infinite_ccw_vertex);

                        Log($"   current_edge supporting line {current_edge.l().l}");
                        Log($"   next_pos_t_ccw {next_pos_t_ccw}");

                        // We can err on the side of having triangulation spokes on the CH boundary,
                        // or we can have collinear points on the CH boundary.
                        // Whichever we pick, we will then have to figure out what the right way
                        // to flip is in one case or another.

                
                        split_this = !( (current_edge).l().l.HasOnPositiveSide(next_pos_t_ccw));
                        if (split_this)
                        {
                            Point2 pos_t_cw;
                            bool infinite_cw_vertex;
                            (pos_t_cw, infinite_cw_vertex) = get_vertex_pos(input, triangle_original_vertex_indices, t_it.t(), cw(t_it.v_in_t_idx()));
                            assert(!infinite_cw_vertex);

                            Log($"   has on pos {current_edge.l().l.HasOnPositiveSide(pos_t_cw)}");
                            Log($"   has on neg {current_edge.l().l.HasOnPositiveSide(pos_t_cw)}");

                            // The triangle created by our splitting the edge from v to ccw(v) will be incident
                            // to the infinite vertex, but the old vertex will not be on the convex hull,
                            // so we need to split things.
                            
                            need_flip = current_edge.l().l.HasOnPositiveSide(pos_t_cw);
                            Log($"   need flip {need_flip}");
                        }
                    }
                }

                if (split_this)
                {
                    Log($" Splitting triangle {t_it.t().Id}");
                    assert(t_it.t().vertex(t_it.v_in_t_idx()) == new_v);
                    new_v_idx++;
                    new_v = new_wavefront_vertices[new_v_idx];

                    Log($"  next vertex is {new_v}");

                    KineticTriangle new_t = split_vertex(t_it.t(), t_it.v_in_t_idx(), current_edge, new_v);
                    assert(current_edge.incident_triangle() == new_t || current_edge_idx  == edges.Count-1);

                    if (need_flip)
                    {
                        Log($" need flip.");
                        Log($"   t  {t_it.t()}");
                        Log($"   tn {new_t}");
                       //  LOG(WARNING) << __FILE__ << ":" << __LINE__ << " " << "untested code path.";
                        new_t.do_raw_flip_inner(0);
                        Log($"   t  {t_it.t()}");
                        Log($"  tn  {new_t}");
                    };

                    t_it = new AroundVertexIterator(new_t, 0);
                    current_edge_idx++;
                    current_edge = edges[current_edge_idx];
                }
                else
                {
                    Log($" Moving on, not splitting this triangle {t_it.t()}");
                    t_it.next();
                }
            }
            assert(new_v_idx == new_wavefront_vertices.Count-1);
            //DBG_FUNC_END(DBG_KT_SETUP);
        }

        ///<summary>
        /// Create bevels
        /// Split vertices for bevels at degree-1 vertices and where we want other reflex vertices to bevel
        /// For each vertex in the triangulation that is not yet initialized, call create_bevels_at_vertex()
        ///</summary>
        private void create_bevels(
        BasicInput input,
        TriangleOriginalVertexIndexList triangle_original_vertex_indices
)
        {
            int initial_triangles_size = triangles.size();

            // verify all existing wavefronts have an incident triangle that
            // reference back to the wavefront.

            foreach (var wf in wavefront_edges)
            {
                assert(wf.incident_triangle() != null);
                assert(wf.incident_triangle().has_wavefront(wf));
            }


            var len = triangles.Count;
            //  set up bevels
            for (var idx =0; idx < len; idx++)
            {
                var t_it = triangles[idx];
                for (int i = 0; i < 3; ++i)
                {
                    if (t_it.vertex(i) != null) continue;
                    assert(t_it.Id < initial_triangles_size);

                    Log($"Doing bevels at t{t_it.Id}; vertex {i}");
                    create_bevels_at_vertex(input, triangle_original_vertex_indices, t_it, i);
                }
            }


            // assert all kinetic triangulation vertices have an assigned
            // wavefront vertex.

            foreach (var t_it in triangles)
            {
                int t_idx = t_it.Id;
                for (int i = 0; i < 3; ++i)
                {
                    assert(t_it.vertex(i) != null);
                }
            }

            // verify all existing wavefronts have an incident triangle that
            // reference back to the wavefront.

            foreach (var wf in wavefront_edges)
            {
                // DBG(DBG_KT_SETUP) << " wf: " << wf;
                assert(wf.incident_triangle() != null);
                assert(wf.incident_triangle().has_wavefront(wf));
            }
        }
        /** return (position, is_infinite) of vertex i of triangle t.
        *
        * This uses the kinetic triangulation's vertex location
        * if already set, and falls back to using the vertex
        * position from the original, underlying, constrained triangulation.
        *
        * The latter only works for triangles that were created initially,
        * not any that are the result of splits (but those should all have
        * kinetic vertices yet).  NO, WRONG, XXX NOT GUARANTEED.
        */

        private partial (Point2, bool) get_vertex_pos(BasicInput input,
          TriangleOriginalVertexIndexList triangle_original_vertex_indices,
          KineticTriangle t,
          int i
        )
        {
            ////DBG_FUNC_BEGIN(//DBG_KT_SETUP);
            //DBG_INDENT_INC();

            WavefrontVertex kv = t.vertex(i);

            Point2 pos = Point2.NaN;
            bool is_inf = false;

            if (kv != null)
            {
                //DBG(//DBG_KT_SETUP) << "  v is " << kv;
                if (kv.is_infinite)
                { /* we already have a vertex set at t's cw */
                    is_inf = true;
                }
                else
                {
                    pos = kv.pos_zero;
                }
            }
            else
            {
                //DBG(//DBG_KT_SETUP) << "  kinetic vertex is not yet set up, using vertex from input triangulation.";
                /* The infinite vertex is associated in kinetic triangles already
                 * in initialize_tds at the very start of setting things up.
                 * So if we are here, it's not the infinite vertex.
                 *
                 * However, it may be a vertex we'll have to bevel later, so it might
                 * still not be set in the kinetic triangulation.  In that case,
                 * get its position from the underlying original triangulation and input.
                 */
                assert(t.Id < triangle_original_vertex_indices.size() / 3); // XXX -- check if this holds
                BasicVertex bv = get_basic_vertex_from_triangle_vertex_indices(input, triangle_original_vertex_indices, t.Id, i);
                pos = bv.p;
            }

            //DBG_INDENT_DEC();
            ////DBG_FUNC_END(//DBG_KT_SETUP);
            return (pos, is_inf);
        }

        private partial void store_initial_wavefront_vertices()
        {
            foreach (var e in wavefront_edges)
            {
                assert(e.is_initial);
                e.set_initial_vertices();
            }
        }
    }
}