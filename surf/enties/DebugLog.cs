#if GODOT
using Godot;
#endif
using System.Diagnostics;

namespace SurfNet

{
    public static class DebugLog
    {

        public static Action<string>? SendWarning = null;
        public static Action<string>? SendMessage = null;
        public static Action<string>? SendRichMessage = null;


        //public static void assert_expensive_eq_ptr(double X, double Y) {
        //        Point_2 __p(X);
        //     Point_2 __q(Y);
        //    assert(Math.Abs(__p.X - __q.X) <= 0 &&   Math.Abs(__p.Y - __q.Y) < 0.0);
        //    }
        

        public static void assert(string message)
        {
            string text = $"{message}";
          
            SendWarning?.Invoke(text);

        }
        public static void assert(bool value, string? message = null)
        {
            if (value) return;

            var st = new StackTrace(1, true); // Saltamos el frame actual (el método IsTrue)
            var frame = st.GetFrame(0);

            string? methodName = frame.GetMethod().Name;
            string? fileName = frame.GetFileName();
            int lineNumber = frame.GetFileLineNumber();

            //Console.ForegroundColor = ConsoleColor.Yellow;
            //Console.WriteLine("⚠️  [CustomAssert Warning]");
            //Console.WriteLine($"Mensaje: {message}");
            //Console.WriteLine($"Ubicación: {fileName}, Método: {methodName}, Línea: {lineNumber}");

            // Mostrar el contenido de la línea que falló (y alrededores)
            if (!string.IsNullOrEmpty(fileName) && File.Exists(fileName))
            {
                string[] lines = File.ReadAllLines(fileName);
                int startLine = Math.Max(0, lineNumber - 3);
                int endLine = Math.Min(lines.Length - 1, lineNumber + 2);

                Console.WriteLine("\n--- Código cerca del fallo ---");
                for (int i = startLine; i <= endLine; i++)
                {
                    string line = lines[i].TrimStart(); // Limpiar espacios iniciales

                    if (i == lineNumber - 1)
                    {
                        //#if GODOT
                        string text = $"{methodName}:l{lineNumber}:{line}:{message}";
                       
                        SendWarning?.Invoke(text);
                    }

                    //#else

                    //       Debug.Fail(line);
                    //#endif

                }
            }
        }


        public static void DotNetDebugAssert(bool condition, string? text = null, string? detail = null)
        {
            Debug.Assert(condition, text, detail);
        }

        public static void DotNetPrintDebugAssert(string text)
        {
            Console.WriteLine(text);
        }


        public static void MetodoEjemplo()
        {
            var stackTrace = new StackTrace();

            for (int i = 0; i < stackTrace.FrameCount; i++)
            {
                var frame = stackTrace.GetFrame(i);
                var method = frame?.GetMethod();

                Console.WriteLine($"Método: {method?.DeclaringType}.{method?.Name}");

                // Si hay símbolos PDB disponibles, podemos obtener archivo y línea
                var fileName = frame?.GetFileName();
                var lineNumber = frame?.GetFileLineNumber();

                if (!string.IsNullOrEmpty(fileName))
                {
                    Console.WriteLine($"  Archivo: {fileName}, Línea: {lineNumber}");
                }
            }
        }
        public static Action<bool, string?, string?> DebugAssertFnc = DotNetDebugAssert;
        public static Action<string> LogPrintFnc = DotNetPrintDebugAssert;

        public static void DebugAssert(bool condition, string? text = null, string? detail = null)
        {
            DebugAssertFnc?.Invoke(condition, text, detail);
        }





        private static bool newLine = true;
        private static int indent = 0;



        public static void LogIndent()
        {
            indent++;

        }
        public static void LogUnindent()
        {
            if (indent > 0)
                indent--;
        }

        public static void Warning(string text = "", ConsoleColor foregroundColor = ConsoleColor.White, ConsoleColor backgroundColor = ConsoleColor.Black)
        {
            Log(text);
        }
        public static void Log(string text = "", ConsoleColor foregroundColor = ConsoleColor.White, ConsoleColor backgroundColor = ConsoleColor.Black)
        {

            SendMessage?.Invoke(text);

            //var bAnterior = Console.BackgroundColor;
            //var fAnterior = Console.ForegroundColor;
            //Console.BackgroundColor = backgroundColor;
            //Console.ForegroundColor = foregroundColor;
            ////  LogPrintFnc?.Invoke($"{new string('\t', indent)}{text}");
            //Console.WriteLine($"{new string(' ', indent * 2)}{text}");
            //Console.BackgroundColor = bAnterior;
            //Console.ForegroundColor = fAnterior;




        }

        public static void LogRich(string text = "")
        {
            SendRichMessage?.Invoke(text);
        }
    }
}
