import bpy
import os
import sys

print("=== Start render_batch.py script ===")

argv = sys.argv
argv = argv[argv.index("--") + 1:] if "--" in argv else []

if len(argv) < 2:
    print("Error: Not enough arguments! Got:", argv)
    sys.exit(1)

obj_path = argv[0]
out_path = argv[1]

print(f"obj_path: {obj_path}")
print(f"out_path: {out_path}")

# 获取water对象
water_obj = bpy.data.objects.get("water")
if water_obj is None:
    print("Error: water object not found in the scene!")
    sys.exit(2)
print("Successfully got water object.")

original_materials = water_obj.data.materials[:]
print(f"Original water materials: {original_materials}")

# 导入新obj mesh
print(f"Importing OBJ mesh from: {obj_path}")
bpy.ops.wm.obj_import(filepath=obj_path)

imported_objs = [obj for obj in bpy.context.selected_objects if obj.type == 'MESH']
if not imported_objs:
    print("Error: No mesh objects found after import!")
    sys.exit(3)

imported_obj = imported_objs[0]
imported_mesh = imported_obj.data
print(f"Imported object: {imported_obj.name}")

# 替换材质
imported_mesh.materials.clear()
for mat in original_materials:
    imported_mesh.materials.append(mat)
print("Materials copied to imported mesh.")

# 替换mesh数据
old_mesh = water_obj.data
water_obj.data = imported_mesh
print("Water object's mesh data has been replaced.")

# 删除导入的临时对象
bpy.data.objects.remove(imported_obj, do_unlink=True)

# 删除无主mesh数据
if old_mesh.users == 0:
    bpy.data.meshes.remove(old_mesh)
    print("Old mesh data removed.")

# 渲染并保存图片
abs_out = os.path.abspath(out_path)
bpy.context.scene.render.filepath = abs_out
print(f"Render output filepath set to: {abs_out}")
print("Starting render...")
bpy.ops.render.render(write_still=True)
print("Render complete. Image saved.")