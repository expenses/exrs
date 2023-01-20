#[macro_use]
extern crate bencher;

extern crate exr;
use exr::prelude::*;

use bencher::Bencher;
use std::io::{Cursor, Write};

fn write_parallel_any_channels_to_buffered(bench: &mut Bencher) {
    let path = "tests/images/valid/custom/crowskull/crow_rle.exr";
    let image = read_all_flat_layers_from_file(path).unwrap();

    bench.iter(||{
        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_parallel_zip1_to_buffered(bench: &mut Bencher) {
    let path = "tests/images/valid/custom/crowskull/crow_rle.exr";

    let mut image = read_first_flat_layer_from_file(path).unwrap();
    image.layer_data.encoding.compression = Compression::ZIP1;

    bench.iter(||{
        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_nonparallel_zip1_to_buffered(bench: &mut Bencher) {
    let path = "tests/images/valid/custom/crowskull/crow_rle.exr";

    let mut image = read_first_flat_layer_from_file(path).unwrap();
    image.layer_data.encoding.compression = Compression::ZIP1;

    bench.iter(||{
        let mut result = Vec::new();
        image.write().non_parallel().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_parallel_zip16_to_buffered(bench: &mut Bencher) {
    let path = "tests/images/valid/custom/crowskull/crow_rle.exr";

    let mut image = read_first_flat_layer_from_file(path).unwrap();
    image.layer_data.encoding.compression = Compression::ZIP16;

    bench.iter(||{
        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_uncompressed_to_buffered(bench: &mut Bencher) {
    let path = "tests/images/valid/custom/crowskull/crow_uncompressed.exr";
    let image = read_all_flat_layers_from_file(path).unwrap();
    assert!(image
        .layer_data
        .iter()
        .all(|layer| layer.encoding.compression == Compression::Uncompressed));

    bench.iter(|| {
        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_scanlines_deinterlaced_custom(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;
    let red = vec![0.2; width * height];
    let green = vec![0.3; width * height];
    let blue = vec![0.6; width * height];

    bench.iter(|| {
        let red = bencher::black_box(red.clone());
        let green = bencher::black_box(green.clone());
        let blue = bencher::black_box(blue.clone());

        let mut result = Vec::new();
        write_exr_fast(
            &mut result,
            width as u32,
            height as u32,
            &red,
            &green,
            &blue,
        );
        bencher::black_box(result);
    })
}

fn write_scanlines_interlaced_custom(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;
    let mut rgba = vec![0.0; width * height * 4];

    for chunk in rgba.chunks_mut(4) {
        chunk[0] = 0.2;
        chunk[1] = 0.3;
        chunk[2] = 0.6;
    }

    bench.iter(|| {
        let rgba = bencher::black_box(rgba.clone());

        let mut red = vec![0.0; width * height];
        let mut green = vec![0.0; width * height];
        let mut blue = vec![0.0; width * height];    

        for (((chunk, red), green), blue) in bencher::black_box(&rgba)
            .chunks_exact(4)
            .zip(&mut red)
            .zip(&mut green)
            .zip(&mut blue)
        {
            *red = chunk[0];
            *green = chunk[1];
            *blue = chunk[2];
        }

        let mut result = Vec::new();
        write_exr_fast(
            &mut result,
            width as u32,
            height as u32,
            bencher::black_box(&red),
            bencher::black_box(&green),
            bencher::black_box(&blue),
        );
        bencher::black_box(result);
    })
}

fn write_scanlines_interlaced_custom_bgr(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;
    let mut rgba = vec![0.0; width * height * 4];

    for chunk in rgba.chunks_mut(4) {
        chunk[0] = 0.2;
        chunk[1] = 0.3;
        chunk[2] = 0.6;
    }

    bench.iter(|| {
        let mut rgba = bencher::black_box(rgba.clone());

        deinterlace_rgba_rows_inplace(&mut rgba, width, height)  ;
        let mut result = Vec::new();
        write_exr_fast_bgr(
            &mut result,
            width as u32,
            height as u32,
            &rgba[..width * height * 3],
        );
        bencher::black_box(result);
    })
}


fn write_scanlines_specificchannels(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;
    let mut rgba = vec![0.0; width * height * 4];

    for chunk in rgba.chunks_mut(4) {
        chunk[0] = 0.2;
        chunk[1] = 0.3;
        chunk[2] = 0.6;
    }

    bench.iter(|| {
        let size = Vec2(width as usize, height as usize);

        let layer1 = Layer::new(
            size,
            LayerAttributes::default(),
            Encoding::UNCOMPRESSED,
            SpecificChannels::rgb(|pos: Vec2<usize>| {
                let offset = (pos.y() * width as usize + pos.x()) * 4;

                (rgba[offset], rgba[offset + 1], rgba[offset + 2])
            }),
        );

        // define the visible area of the canvas
        let attributes = ImageAttributes::new(
            // the pixel section that should be shown
            IntegerBounds::from_dimensions(size),
        );

        let image = Image::empty(attributes).with_layer(layer1); // add an rgba layer of different type, `SpecificChannels<ClosureB>`, not possible with a vector

        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_scanlines_deinterlaced_anychannels(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;

    let size = Vec2(width as usize, height as usize);

    let red = vec![0.2; width * height];
    let green = vec![0.3; width * height];
    let blue = vec![0.6; width * height];

    bench.iter(|| {
        let red = bencher::black_box(red.clone());
        let green = bencher::black_box(green.clone());
        let blue = bencher::black_box(blue.clone());

        let red_s = exr::image::FlatSamples::F32(red);
        let green_s = exr::image::FlatSamples::F32(green);
        let blue_s = exr::image::FlatSamples::F32(blue);

        // define the visible area of the canvas
        let attributes = ImageAttributes::new(
            // the pixel section that should be shown
            IntegerBounds::from_dimensions(size),
        );

        let layer1 = Layer::new(
            size,
            LayerAttributes::default(),
            Encoding::UNCOMPRESSED,
            exr::image::AnyChannels::sort(smallvec::smallvec![
                exr::image::AnyChannel::new("R", red_s),
                exr::image::AnyChannel::new("G", green_s),
                exr::image::AnyChannel::new("B", blue_s),
            ]),
        );

        let image = Image::empty(attributes).with_layer(layer1); // add an rgba layer of different type, `SpecificChannels<ClosureB>`, not possible with a vector

        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_scanlines_interlaced_anychannels(bench: &mut Bencher) {
    let width = 3000;
    let height = 3000;

    let size = Vec2(width as usize, height as usize);

    let mut rgba = vec![0.0_f32; width * height * 4];

    for chunk in rgba.chunks_mut(4) {
        chunk[0] = 0.2;
        chunk[1] = 0.3;
        chunk[2] = 0.6;
    }

    bench.iter(|| {
        let rgba = bencher::black_box(rgba.clone());

        let mut red = vec![0.0; width * height];
        let mut green = vec![0.0; width * height];
        let mut blue = vec![0.0; width * height];    

        for (((chunk, red), green), blue) in bencher::black_box(&rgba)
            .chunks_exact(4)
            .zip(&mut red)
            .zip(&mut green)
            .zip(&mut blue)
        {
            *red = chunk[0];
            *green = chunk[1];
            *blue = chunk[2];
        }

        let red_s = exr::image::FlatSamples::F32(red);
        let green_s = exr::image::FlatSamples::F32(green);
        let blue_s = exr::image::FlatSamples::F32(blue);

        // define the visible area of the canvas
        let attributes = ImageAttributes::new(
            // the pixel section that should be shown
            IntegerBounds::from_dimensions(size),
        );

        let layer1 = Layer::new(
            size,
            LayerAttributes::default(),
            Encoding::UNCOMPRESSED,
            exr::image::AnyChannels::sort(smallvec::smallvec![
                exr::image::AnyChannel::new("R", &red_s),
                exr::image::AnyChannel::new("G", &green_s),
                exr::image::AnyChannel::new("B", &blue_s),
            ]),
        );

        let image = Image::empty(attributes).with_layer(layer1); // add an rgba layer of different type, `SpecificChannels<ClosureB>`, not possible with a vector

        let mut result = Vec::new();
        image.write().to_buffered(Cursor::new(&mut result)).unwrap();
        bencher::black_box(result);
    })
}

fn write_exr_fast<W: Write>(
    output: &mut W,
    width: u32,
    height: u32,
    red_floats: &[f32],
    green_floats: &[f32],
    blue_floats: &[f32],
) {
    let mut output_vec = Vec::new();

    let size = exr::math::Vec2(width as usize, height as usize);

    let num_channels = 3;
    let subpixel_width = 4;
    let pixel_width = num_channels * subpixel_width;

    let header = exr::meta::header::Header {
        channels: exr::meta::attribute::ChannelList::new(
            smallvec::smallvec![
                exr::meta::attribute::ChannelDescription::named(
                    "B",
                    exr::meta::attribute::SampleType::F32
                ),
                exr::meta::attribute::ChannelDescription::named(
                    "G",
                    exr::meta::attribute::SampleType::F32
                ),
                exr::meta::attribute::ChannelDescription::named(
                    "R",
                    exr::meta::attribute::SampleType::F32
                ),
            ]
            .into(),
        ),
        compression: exr::compression::Compression::Uncompressed,
        blocks: exr::meta::BlockDescription::ScanLines,
        line_order: exr::meta::attribute::LineOrder::Increasing,
        layer_size: size,
        deep: false,
        deep_data_version: None,
        chunk_count: height as usize,
        max_samples_per_pixel: None,
        shared_attributes: exr::meta::header::ImageAttributes::with_size(size),
        own_attributes: Default::default(),
    };

    exr::meta::MetaData::write_validating_to_buffered(&mut output_vec, &[header], true).unwrap();

    output.write_all(&output_vec).unwrap();

    let mut offset_tables = Vec::with_capacity(height as usize);

    let base_offset = output_vec.len() as u64 + height as u64 * 8;

    drop(output_vec);

    for i in 0..height as u64 {
        offset_tables.push(base_offset as u64 + i * ((width * pixel_width) + 8) as u64);
    }

    output.write_all(bytemuck::cast_slice(&offset_tables));

    for y in 0..height {
        output.write_all(&y.to_le_bytes());
        output.write_all(&(width * pixel_width).to_le_bytes());

        let y = y as usize;
        let width = width as usize;

        output.write_all(bytemuck::cast_slice(
            &blue_floats[y * width..(y + 1) * width],
        ));
        output.write_all(bytemuck::cast_slice(
            &green_floats[y * width..(y + 1) * width],
        ));
        output.write_all(bytemuck::cast_slice(
            &red_floats[y * width..(y + 1) * width],
        ));
    }
}


fn write_exr_fast_bgr<W: Write>(
    output: &mut W,
    width: u32,
    height: u32,
    bgr: &[f32],
) {
    let mut output_vec = Vec::new();

    let size = exr::math::Vec2(width as usize, height as usize);

    let num_channels = 3;
    let subpixel_width = 4;
    let pixel_width = num_channels * subpixel_width;

    let header = exr::meta::header::Header {
        channels: exr::meta::attribute::ChannelList::new(
            smallvec::smallvec![
                exr::meta::attribute::ChannelDescription::named(
                    "B",
                    exr::meta::attribute::SampleType::F32
                ),
                exr::meta::attribute::ChannelDescription::named(
                    "G",
                    exr::meta::attribute::SampleType::F32
                ),
                exr::meta::attribute::ChannelDescription::named(
                    "R",
                    exr::meta::attribute::SampleType::F32
                ),
            ]
            .into(),
        ),
        compression: exr::compression::Compression::Uncompressed,
        blocks: exr::meta::BlockDescription::ScanLines,
        line_order: exr::meta::attribute::LineOrder::Increasing,
        layer_size: size,
        deep: false,
        deep_data_version: None,
        chunk_count: height as usize,
        max_samples_per_pixel: None,
        shared_attributes: exr::meta::header::ImageAttributes::with_size(size),
        own_attributes: Default::default(),
    };

    exr::meta::MetaData::write_validating_to_buffered(&mut output_vec, &[header], true).unwrap();

    output.write_all(&output_vec).unwrap();

    let mut offset_tables = Vec::with_capacity(height as usize);

    let base_offset = output_vec.len() as u64 + height as u64 * 8;

    drop(output_vec);

    for i in 0..height as u64 {
        offset_tables.push(base_offset as u64 + i * ((width * pixel_width) + 8) as u64);
    }

    output.write_all(bytemuck::cast_slice(&offset_tables));

    for y in 0..height {
        output.write_all(&y.to_le_bytes());
        output.write_all(&(width * pixel_width).to_le_bytes());

        let y = y as usize;
        let width = width as usize;

        output.write_all(bytemuck::cast_slice(
            &bgr[y * width * 3..(y + 1) * width * 3],
        ));
    }
}


benchmark_group!(write,
    write_parallel_any_channels_to_buffered,
    write_nonparallel_zip1_to_buffered,
    write_parallel_zip1_to_buffered,
    write_parallel_zip16_to_buffered,
    write_uncompressed_to_buffered,
    write_scanlines_deinterlaced_custom,
    write_scanlines_interlaced_custom,
    write_scanlines_interlaced_anychannels,
    write_scanlines_deinterlaced_anychannels,
    write_scanlines_specificchannels,
    write_scanlines_interlaced_custom_bgr,
);

benchmark_main!(write);

fn deinterlace_rgba_rows_inplace(rgba: &mut [f32], width: usize, height: usize) {
    let mut row = vec![0.0; width * 3];

    for y in 0 .. height {
        let mut i = 0;
        {
            let input = &rgba[y * width * 4 .. (y + 1) * width * 4];
            let (mut blue, mut remaining) = row.split_at_mut(width);
            let (mut green, mut red) = remaining.split_at_mut(width);
            for (((chunk, red), green), blue) in input.chunks_exact(4).zip(red).zip(green).zip(blue) {
                *red = chunk[0];
                *green = chunk[1];
                *blue = chunk[2];
            }
        }
        rgba[y * width * 3..(y +1) * width * 3].copy_from_slice(&row);
    }
}
