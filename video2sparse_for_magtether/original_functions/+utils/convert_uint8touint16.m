function data = convert_uint8touint16(data)
    data = uint16(double(data)./(2^8-1).*(2^16-1));
end