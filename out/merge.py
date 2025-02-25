from PIL import Image

def merge_images(image_paths, direction='horizontal'):
    """
    合并多张图片为一张图片。
    
    参数:
    image_paths: 图片路径的列表
    direction: 合并方向，'horizontal' 或 'vertical'
    
    返回:
    合并后的图片
    """
    if direction == 'horizontal':
        # 计算合并后图片的宽度和高度
        total_width = sum(Image.open(img).size[0] for img in image_paths)
        max_height = max(Image.open(img).size[1] for img in image_paths)
    elif direction == 'vertical':
        total_height = sum(Image.open(img).size[1] for img in image_paths)
        max_width = max(Image.open(img).size[0] for img in image_paths)
    else:
        raise ValueError("方向只能是 'horizontal' 或 'vertical'")

    # 创建一个新的空白图片
    to_image = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    y_offset = 0

    for img_path in image_paths:
        with Image.open(img_path) as img:
            if direction == 'horizontal':
                to_image.paste(img, (x_offset, 0))
                x_offset += img.size[0]
            elif direction == 'vertical':
                to_image.paste(img, (0, y_offset))
                y_offset += img.size[1]

    return to_image

out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'

methods = ['regular', 'whitenoise', 'dartthrowing', 'stratified', 'poissondisk', 'NESOTS']

image_paths = [out_dir + f'show_{m}.jpg' for m in methods]

merged_image = merge_images(image_paths, direction='horizontal')
merged_image.save(out_dir+'merged_image.jpg')