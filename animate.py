from utilities import make_mp4, make_gif


stations = [10,20,30,40,50,60,70,80,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540]

# make_mp4(
#       [f"comparison_x={x}.png" for x in stations],
#       "du.mp4",
#       fps=7,
#       quality=10,
#   )

make_gif(
      [f"S800_dzTauxz_x={x}.png" for x in stations],
      "dztauxz-full.gif",
      fps=7,
  )

# make_gif(
#       [f"S800_x={x}.png" for x in stations],
#       "S800.gif",
#       fps=7,
#   )