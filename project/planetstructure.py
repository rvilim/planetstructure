import phases
import params
import structure


def main():
    recipe_file = "recipes/earth.json"

    params.get_params('params/params.ini', recipe_file)

    structure.structure(6371000.0)

if __name__ == "__main__":
    main()
