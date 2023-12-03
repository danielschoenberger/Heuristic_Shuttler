mixed_list = [1, (2, 3), (4, 2), (5, 6), 7]
value_to_check = 4

found = any(
    value_to_check in element if isinstance(element, tuple) else element == value_to_check for element in mixed_list
)

if found:
    print(f"{value_to_check} is present in the mixed list.")
else:
    print(f"{value_to_check} is not present in the mixed list.")
