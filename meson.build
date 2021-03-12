project('tiny_ising', 'c', default_options: [ 'warning_level=2' ])

cc = meson.get_compiler('c')
m_dep = cc.find_library('m')
gomp_dep = cc.find_library('gomp')

tiny_ising = executable('tiny_ising',
    [ 'tiny_ising.c', 'ising.c' ],
    dependencies: [ m_dep, gomp_dep ]
)

gl_dep = cc.find_library('GL', required: false)
glfw_dep = cc.find_library('glfw', required: false)

demo = executable('demo',
    [ 'demo.c', 'ising.c' ],
    dependencies: [ m_dep, gomp_dep, glfw_dep, gl_dep ]
)
