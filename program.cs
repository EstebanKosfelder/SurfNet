﻿using SurfNet;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;



SkeletonStructure skeletonStructure = new SkeletonStructure(null);

 var builder =  skeletonStructure.kt;



var polys = KineticTriangulation.CargarListaVector2("F:\\trabajos\\g_test.json");

foreach (var poly in polys)
{

    builder.EnterContour(poly.Select(p => new Point2(p.X, p.Y)));
}



builder.Initialize();



var wp = skeletonStructure.wp;

wp.setup_queue(builder);

while (!wp.propagation_complete() && !builder.queue.empty())
{
    //var item =(EventQueueItem) ( builder.EventQueue.pop();
    //builder.handle_event(item.priority);




    wp.advance_step();

}



Console.WriteLine("fin");
