from PIL import Image



out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'

methods = ['regular', 'whitenoise', 'dartthrowing', 'stratified', 'poissondisk', 'NESOTS', 'HCPSOTS']

image_paths = [out_dir + f'show_{m}.jpg' for m in methods]

# 获取图片尺寸（这里假设所有图片尺寸相同，如果不同需要调整）
image = Image.open(image_paths[0])
image_width, image_height = image.size

# 设置画布尺寸（根据图片数量和布局设置）
canvas_width = 4 * image_width  # 假设每张图片宽度相同，4张一行
canvas_height = 2 * image_height  # 两行
# 创建空白画布
canvas = Image.new('RGB', (canvas_width, canvas_height), 'white')

# 粘贴图片到画布上
x_offset = 0
y_offset = 0
for i, image_path in enumerate(image_paths):
    image = Image.open(image_path)
    if i < 4:
        # 第一行
        canvas.paste(image, (x_offset, y_offset))
        x_offset += image_width
    else:
        # 重置x偏移量，准备粘贴到第二行
        if x_offset == canvas_width:
            x_offset = int(0.5*image_width)
            y_offset += image_height
        # 第二行
        canvas.paste(image, (x_offset, y_offset))
        x_offset += image_width

canvas.save(out_dir+'merged_sphere_noise.jpg')

