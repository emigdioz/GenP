{
  "targets": [
    {
      "target_name": "nodeGP",
      "sources": [ "../fire.cpp", "../exchange.cpp", "addon.cpp"],
      "cflags": ["-Wall", "-std=c++11","-fexceptions"],
      "cflags_cc": [ '-fexceptions' ],
      "libraries" : ["-lGP-parallel -lOpenCL -lopenblas"],
      "include_dirs" : ["../","<!(nodejs -e \"require('nan')\")"],
      "conditions": [
        [ 'OS=="mac"', {
            "xcode_settings": {
                'OTHER_CPLUSPLUSFLAGS' : ['-std=c++11','-stdlib=libc++'],
                'OTHER_LDFLAGS': ['-stdlib=libc++'],
                'MACOSX_DEPLOYMENT_TARGET': '10.7' }
            }
        ]
      ]
    }
  ]
}
