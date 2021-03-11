project('tiny_ising', 'c')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m')
gomp_dep = cc.find_library('gomp')

tiny_ising = executable('tiny_ising',
    'tiny_ising.c',
    dependencies: [ m_dep, gomp_dep ]
)

gl_dep = cc.find_library('GL', required: false)
glfw_dep = cc.find_library('glfw', required: false)

demo = executable('demo',
    'tiny_ising.c',
    c_args: '-DGRAPHICS',
    dependencies: [ m_dep, gomp_dep, glfw_dep, gl_dep ]
)
