project('tiny_ising', 'c',
        default_options: [ 'c_std=c11', 'warning_level=2' ]
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m')

tiny_ising = executable('tiny_ising',
    [ 'tiny_ising.c', 'ising.c', 'wtime.c' ],
    dependencies: [ m_dep ]
)

gl_dep = cc.find_library('GL', required: false)
glfw_dep = cc.find_library('glfw', required: false)

demo = executable('demo',
    [ 'demo.c', 'ising.c', 'wtime.c' ],
    dependencies: [ m_dep, glfw_dep, gl_dep ]
)
