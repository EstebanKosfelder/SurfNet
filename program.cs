using SurfNet;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;



SkeletonStructure skeletonStructure = new SkeletonStructure(null);

 var builder =  skeletonStructure.kt;


DebugLog.SendWarning = (s) =>
{
    var a = Console.ForegroundColor;
    Console.ForegroundColor = ConsoleColor.Yellow;
    Console.WriteLine(s);
    Console.ForegroundColor = a;
};

DebugLog.SendMessage = (s) =>
{
    Console.WriteLine(s);
 
};

 var polys = KineticTriangulation.CargarListaVector2("D:\\Documents\\Trabajos\\test1.json");
//var polys = KineticTriangulation.CargarListaVector2("D:\\Documents\\Trabajos\\g2.json");

foreach (var poly in polys)
{

    builder.EnterContour(poly.Select(p => new Point2(Math.Round(p.X,3), Math.Round(p.Y))));
}




skeletonStructure.initialize();



var wp = skeletonStructure.wp;



while (!wp.propagation_complete() && !builder.queue.empty())
{
    //var item =(EventQueueItem) ( builder.EventQueue.pop();
    //builder.handle_event(item.priority);




    wp.advance_step();

}



Console.WriteLine("fin");
